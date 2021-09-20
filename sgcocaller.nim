
## two modules (generate pair-wise fragments, and then run viterbi)


import os
import docopt
import strutils
import hts
import tables
import sequtils
import private/utils
import private/fragments
import math
import streams
import private/graph
import private/findpath
import private/blocks_utils
import private/readfmf
import private/phase_blocks


proc getfmf(ibam:Bam, ivcf:VCF, barcodeTable:TableRef, outfmf:FileStream,maxTotalReads:int,
                  minTotalReads:int,  mapq: int,minbsq:int,mindp:int,barcodeTag:string, minsnpdepth:int):int = 
  var variantIndex = 1
  var cellFragmentsTable = newTable[string,Fragment]()

  echo "generated fragments from " & $barcodeTable.len & " single cells"

  for rec in ivcf.query("*"): 
    var cellGtNodes = findGtNodes(rec=rec, variantIndex= variantIndex, ibam=ibam, maxTotalReads=maxTotalReads,
                                  minTotalReads=minTotalReads, mapq = mapq, barcodeTable = barcodeTable, minbsq=minbsq,barcodeTag=barcodeTag,
                                  minCellDp = mindp, minSnpDepth = minsnpdepth)
    cellFragmentsTable = updateCellFragment(barcodedNodesTable=cellGtNodes,cellFragmentsTable = cellFragmentsTable, fmf_out = outfmf)              
    variantIndex = variantIndex + 1
  echo "processed " & $(variantIndex) & " variants"
  return 0
proc phaseBlocks(fmf_file:string, inputVCF:string, outputVCF:string, block_hapfile:string, threads:int):int =
  var cellGenoTable: TableRef[string,TableRef[string,int]]
  var cellBlockLeftPhaseTable,cellBlockRightPhaseTable :TableRef[string, seq[int]]
  var chr_blocks:seq[Block]
  var blockkseq:seq[Block_linkage]
  var phased_blocks: seq[int]
  echo "Phasing for hapfile " & block_hapfile
  cellGenoTable = newTable[string,TableRef[string,int]]()
  discard readFmf(fmf_file,cellGenoTable)
  chr_blocks = read_blocks(block_hapfile)
  cellBlockLeftPhaseTable = newTable[string, seq[int]]()
  cellBlockRightPhaseTable = newTable[string, seq[int]]()
  ## blocks at the ends with 20 SNPs or fewer, removed?
  ## blocks not at the ends with 20SNPs or fewer are not used for linking blocks but if their linkage with pre block and post-block is consistent, 
  ## these short blocks' haplotypes are retained. 
  ## 
  ## What if blocks cannot be linked ? the number of 11(00) == 10(01)
  echo "Total blocks " & $chr_blocks.len
  var contig_blocks:seq[int]
  for b_j, bl in chr_blocks:
    echo "block " & $b_j & " len: " & $bl.h0.len
    # echo "last 5" & $bl.h0[(high(bl.h0)-5)..high(bl.h0)]
    # echo "first 5" & $bl.h0[0..5]
    discard block_phase_in_cells(cellGenoTable, bl.snpPos, bl.h0, true, cellBlockLeftPhaseTable)
    discard block_phase_in_cells(cellGenoTable, bl.snpPos, bl.h0, false, cellBlockRightPhaseTable)
    contig_blocks = find_contig_blocks(chr_blocks, cellBlockLeftPhaseTable, cellBlockRightPhaseTable)
  # echo "blocks for contigs" & $contig_blocks
  blockkseq = build_contig(contig_blocks = contig_blocks, cellBlockLeftPhaseTable, cellBlockRightPhaseTable)
  for bk in blockkseq:
    echo "prevBlock" &  $bk.prevBlock & "nextBlock" &  $bk.nextBlock & " 00 type " & $ $bk.linkageType_00 &  " 01 type " & $bk.linkageType_01 &  " switch posterior prob " & $bk.switch_posterior_prob  & " switch value :" & $bk.switch
  phased_blocks = infer_non_contig_blocks(contig_blocks, blockkseq, chr_blocks, cellBlockLeftPhaseTable, cellBlockRightPhaseTable)

  # echo phased_blocks

  discard  write_linked_blocks(inputVCF,outputVCF, blocks = chr_blocks, blocks_phase = phased_blocks, threads= threads)
  return 0
proc sgcocaller(threads:int, ivcf:VCF, barcodeTable:TableRef,
                ibam:Bam, out_dir:string, mapq:int, 
                minbsq:int, mintotal:int, maxtotal:int, mindp:int, maxdp:int, 
                thetaREF:float, thetaALT:float, cmPmb:float,s_Chrs:seq,barcodeTag:string,phased:bool): int =
  
  var bulkBam = false
  if barcodeTable.len == 1 and barcodeTable.hasKey("bulk"):
    echo "running in bulk mode and all reads are regarded as from one sample/cell"
    bulkBam = true
  else:
    echo "running sgcoaller for " & $barcodeTable.len & " single cells"

  let initProb:array[stateRef..stateAlt, float]=[0.5,0.5]
  ## iterate through each selected chromosome
  for chrom in s_Chrs:
    var snpIndex, nnsize = 0
    ## number of non zeros
    var scSpermSeq:SeqSpermViNodes
    ## matches with the order in barcodeTable
    scSpermSeq.setLen(barcodeTable.len)
    var outFileSNPanno,outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx,viSegmentInfo:FileStream

    try:
      outFileSNPanno = openFileStream(out_dir & chrom & "_snpAnnot.txt", fmWrite)
      outFileTotalCountMtx = openFileStream(out_dir & chrom & "_totalCount.mtx", fmWrite)
      outFileAltCountMtx = openFileStream(out_dir & chrom & "_altCount.mtx", fmWrite)
      outFileVStateMtx = openFileStream(out_dir & chrom & "_vi.mtx", fmWrite)
      viSegmentInfo = openFileStream(out_dir & chrom & "_viSegInfo.txt", fmWrite)
    except:
      stderr.write getCurrentExceptionMsg()
    
    ## write headers to those mtx files and the first line place holder for total_row total_column total_entry
    let sparseMatrixHeader = "%%MatrixMarket matrix coordinate integer general"
    
    for outFileStream in [outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx]:
      outFileStream.writeLine(sparseMatrixHeader)
      outFileStream.writeLine(' '.repeat(50))
    
    outFileSNPanno.writeLine(join(["POS","REF", "ALT"], sep = "\t"))
    var ac = new_seq[int32](10)
    for rec in ivcf.query(chrom):
      if rec.ALT.len > 1 or rec.ALT.len == 0 : continue
      ## alleleCountTable contains for this rec.POS, each cell barcode's allele counts 
      var rec_alt = rec.ALT[0][0]
      var rec_ref = rec.REF[0]
      if phased:
        var f = rec.format()
        var gts = f.genotypes(ac)
        if ac[0] == 4 and ac[1] == 3:
          # gt is 1|0
            rec_alt = rec.REF[0]
            rec_ref = rec.ALT[0][0]
          # if not gt == 0|1 continue
        elif not (ac[0] == 2 and ac[1] == 5): continue
      var alleleCountTable = countAllele(ibam=ibam, chrom=chrom, mapq=mapq, barcodeTable=barcodeTable,minbsq=minbsq,
                                        maxTotalReads = maxtotal, minTotalReads = mintotal,bulkBam = bulkBam, barcodeTag = barcodeTag, startPos = rec.POS.cint-1, stopPos=rec.POS.cint,rec_alt =  rec_alt, rec_ref = rec_ref)
      if alleleCountTable.len==0: continue

      ## add to snpAnnoSeq, later write to SNPannot file, which contains SNP.pos, SNP.ref,SNP.alt; The rowAnnotations
      snpIndex += 1
      outFileSNPanno.writeLine(join([$rec.POS, $rec_ref, $rec_alt], sep="\t") )
      discard addViNode(barcodeTable = barcodeTable,  alleleCountTable = alleleCountTable,  scSpermSeq = scSpermSeq,
               outFileTotalCountMtx = outFileTotalCountMtx, outFileAltCountMtx  = outFileAltCountMtx, nnsize = nnsize,
               mindp = mindp, maxdp = maxdp, thetaRef = thetaRef, thetaAlt = thetaAlt, snpIndex = snpIndex, rec_pos=int(rec.POS),
               initProb = initProb, cmPmb = cmPmb)
    var posEnd: int64
    var inferProb,reverseProb = 0.0
    var currentEm: array[stateRef..stateAlt, float]
    var lastNode: ViNode
    var spermVNseq: SpermViNodes
    var ithSNP: int

    for ithSperm in 0..(scSpermSeq.len-1):
      ## rightGap,leftGap = [0.0,0.0]
      spermVNseq = scSpermSeq[ithSperm]
      if spermVNseq.viNodeseq.len==0: continue
      lastNode = spermVNseq.viNodeseq[high(spermVNseq.viNodeseq)]
      currentEm = getEmission(thetaRef=thetaRef,thetaAlt=thetaAlt,cRef=lastNode.cRef, cAlt=lastNode.cAlt)
      if lastNode.pathScore[stateRef] > lastNode.pathScore[stateAlt]:
        scSpermSeq[ithSperm].viNodeseq[high(scSpermSeq[ithSperm].viNodeseq)].state = stateRef
        ithSNP = scSpermSeq[ithSperm].spermSnpIndexLookUp[high(scSpermSeq[ithSperm].viNodeseq)+1]
        outFileVStateMtx.writeLine($ithSNP & " " & $(ithSperm+1) & " 1")
        inferProb  = currentEm[stateRef]
        reverseProb = currentEm[stateAlt]
      else:
        scSpermSeq[ithSperm].viNodeseq[high(scSpermSeq[ithSperm].viNodeseq)].state = stateAlt
        ithSNP = scSpermSeq[ithSperm].spermSnpIndexLookUp[high(scSpermSeq[ithSperm].viNodeseq)+1]
        outFileVStateMtx.writeLine($ithSNP & " " & $(ithSperm+1) & " 2")
        inferProb  = currentEm[stateAlt]
        reverseProb = currentEm[stateRef]
      posEnd = lastNode.pos
      ## call pathTrackBack will write the most probably hidden state seq and viterbi segment info to the  relevant File Streams.
      discard pathTrackBack(currentSperm = scSpermSeq[ithSperm], ithSperm = ithSperm,  thetaRef = thetaRef, 
                           thetaAlt = thetaAlt, cmPmb = cmPmb,outFileVStateMtx = outFileVStateMtx,
                           viSegmentInfo = viSegmentInfo, posEnd = posEnd,  inferProb = inferProb,reverseProb = reverseProb)
    for outFileStream in [outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx]:
      outFileStream.setPosition(49)
      outFileStream.write($snpIndex & " " & $barcodeTable.len & " " & $nnsize) 
    for outFileStream in [outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx,outFileSNPanno,viSegmentInfo]:
      outFileStream.close()
  return 0

when(isMainModule):
  let version = "0.3.3"
  var doc = format("""
  $version

  Usage:
      sgcocaller fmf [options] <BAM> <VCF> <barcodeFile> <out_prefix>
      sgcocaller xo [options] <BAM> <VCF> <barcodeFile> <out_prefix>
      sgcocaller phase_blocks [options] <VCF> <fmf> <hapfile> <out_vcf>
      

Arguments:

  <BAM> the read alignment file with records of single-cell DNA reads
  
  <VCF> the variant call file with records of SNPs

  <barcodeFile> the text file containing the list of cell barcodes

  <out_prefix>  the prefix of output files

  <fmf> the fragment file from running sgcocaller fmf

  <out_vcf> the output vcf aftering phasing blocks in hapfile


Options:
  -t --threads <threads>  number of BAM decompression threads [default: 4]
  --barcodeTag <barcodeTag>  the cell barcode tag in BAM [default: CB]
  --minMAPQ <mapq>  Minimum MAPQ for read filtering [default: 20]
  --baseq <baseq>  base quality threshold for a base to be used for counting [default: 13]
  --chrom <chrom>  the selected chromsome (whole genome if not supplied,separate by comma if multiple chroms)
  --minDP <minDP>  the minimum DP for a SNP to be included in the output file [default: 1]
  --maxDP <maxDP>  the maximum DP for a SNP to be included in the output file [default: 5]
  --maxTotalDP <maxTotalDP>  the maximum DP across all barcodes for a SNP to be included in the output file [default: 25]
  --minTotalDP <minTotalDP>  the minimum DP across all barcodes for a SNP to be included in the output file [default: 10]
  --minSNPdepth <minSNPdepth>  the minimum depth of coverage for a SNPs to be includes in generated fragments [default: 2]
  --thetaREF <thetaREF>  the theta for the binomial distribution conditioning on hidden state being REF [default: 0.1]
  --thetaALT <thetaALT>  the theta for the binomial distribution conditioning on hidden state being ALT [default: 0.9]
  --cmPmb <cmPmb>  the average centiMorgan distances per megabases default 0.1 cm per Mb [default: 0.1]
  --phased  the VCF contains the phased GT of heterozygous
  -h --help  show help


  Examples
      ./sgcocaller fmf --threads 4 AAAGTAGCACGTCTCT-1.raw.bam AAAGTAGCACGTCTCT-1.raw.bam.dp3.alt.vcf.gz barcodeFile.tsv ./percell/cellfragment.fmf
      ./sgcocaller xo --threads 4 AAAGTAGCACGTCTCT-1.raw.bam AAAGTAGCACGTCTCT-1.raw.bam.dp3.alt.vcf.gz barcodeFile.tsv ./percell/ccsnp
      ./sgcocaller phase_blocks --threads 4 AAAGTAGCACGTCTCT-1.raw.bam.dp3.alt.vcf.gz ./percell/cellfragment.fmf ./percell/ccsnp/linkedBlocks.vcf.gz

""" % ["version", version])

  let args = docopt(doc, version=version)
  
  var
    threads:int
    selectedChrs:string
    out_dir:string
    mapq:int
    minbsq:int
    mintotal:int
    maxtotal:int
    mindp:int
    maxdp:int
    thetaREF:float
    thetaALT:float
    cmPmb:float
    barcodeTag="CB"
    minsnpdepth:int
    out_vcf,fmf,barcodeFile,bamfile,vcff,hapfile:string
    vcfGtPhased = false

  echo $args
  threads = parse_int($args["--threads"])
  barcodeTag = $args["--barcodeTag"]
  mindp = parse_int($args["--minDP"])
  maxdp = parse_int($args["--maxDP"])
  maxtotal = parse_int($args["--maxTotalDP"])
  mintotal = parse_int($args["--minTotalDP"])
  minsnpdepth = parse_int($args["--minSNPdepth"])
  mapq = parse_int($args["--minMAPQ"])
  minbsq = parse_int($args["--baseq"])
  thetaRef = parse_float($args["--thetaREF"])
  thetaAlt = parse_float($args["--thetaALT"])
  cmPmb = parse_float($args["--cmPmb"])
  vcfGtPhased = parse_bool($args["--phased"])
  
  if args["fmf"] or args["xo"]:
    var 
      ibam:Bam
      ivcf:VCF
      s_Chrs: seq[string]
    bamfile = $args["<BAM>"]
    barcodeFile = $args["<barcodeFile>"]
    out_dir = $args["<out_prefix>"]
    vcff = $args["<VCF>"]

    if not open(ibam, bamfile, threads=threads, index = true):
        quit "couldn't open input bam"
    if not open(ivcf, vcff, threads=threads):
        quit "couldn't open: vcf file"

    if($args["--chrom"] != "nil"):
      selectedChrs = $args["--chrom"]
      s_Chrs = selectedChrs.split(',')
    else:
      ## find all chroms from headers of vcf file (Contigs)
      let contigs = ivcf.contigs 
      s_Chrs = map(contigs, proc(x: Contig): string = x.name)

    var hf = hts.hts_open(cstring(barcodeFile), "r")
    #### TODO : Table size
    var barcodeTable =  newTable[string,int](initialSize = 1024)
    var kstr: hts.kstring_t
    kstr.l = 0
    kstr.m = 0
    kstr.s = nil
    var ithSperm = 0
    ## initiate the table with CB as keys, allele counts (ref object) as elements
    while hts_getline(hf, cint(10), addr kstr) > 0:
      if $kstr.s[0] == "#":
        continue
      var v = $kstr.s
      discard barcodeTable.hasKeyOrPut(v, ithSperm)
      ithSperm.inc
    discard hf.hts_close()
    if args["fmf"]:
      echo "generating fragment file to " & out_dir
      var outfmf:FileStream
      try: 
        outfmf = openFileStream(out_dir, fmWrite)
      except:
        stderr.write getCurrentExceptionMsg()
      discard getfmf(ibam = ibam, ivcf = ivcf, barcodeTable = barcodeTable, outfmf = outfmf, maxTotalReads = maxtotal,
                    minTotalReads = mintotal, mapq = mapq, minbsq =minbsq, mindp = mindp, barcodeTag = barcodeTag,minsnpdepth=minsnpdepth)
      close(outfmf)
    if args["xo"]:
      echo "running crossover calling \n"
      echo "running for these chromosomes as specified in the header of the provided VCF file " & $s_Chrs
      discard sgcocaller(threads, ivcf, barcodeTable, ibam,
                      out_dir, mapq, minbsq, mintotal, 
                      maxtotal, mindp, maxdp, thetaREF, thetaALT, cmPmb,s_Chrs,barcodeTag,phased = vcfGtPhased)
    ibam.close()
    ivcf.close()
  if args["phase_blocks"]:
    fmf = $args["<fmf>"]
    out_vcf = $args["<out_vcf>"]
    vcff = $args["<VCF>"]
    hapfile =  $args["<hapfile>"]
    echo "running phasing blocks \n"
    discard phaseBlocks(fmf_file = fmf, inputVCF=vcff, outputVCF=out_vcf, block_hapfile=hapfile, threads=threads)


 