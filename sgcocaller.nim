
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

proc getfmf(ibam:Bam, ivcf:VCF, barcodeTable:OrderedTableRef, outfmf:FileStream,maxTotalReads:int,
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

proc sgcocaller(threads:int, ivcf:VCF, barcodeTable:OrderedTableRef,
                ibam:Bam, out_dir:string, mapq:int, 
                minbsq:int, mintotal:int, maxtotal:int, mindp:int, maxdp:int, 
                thetaREF:float, thetaALT:float, cmPmb:float,s_Chrs:seq,barcodeTag:string): int =
  
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

    for rec in ivcf.query(chrom):
      if rec.ALT.len > 1 or rec.ALT.len == 0 : continue
      ## alleleCountTable contains for this rec.POS, each cell barcode's allele counts 
      var alleleCountTable = countAllele(rec=rec, ibam=ibam, chrom=chrom, mapq=mapq, barcodeTable=barcodeTable,minbsq=minbsq,
                                        maxTotalReads = maxtotal, minTotalReads = mintotal,bulkBam = bulkBam, barcodeTag = barcodeTag)
      if alleleCountTable.len==0: continue
      var rec_alt:char
      rec_alt = rec.ALT[0][0]
      ## add to snpAnnoSeq, later write to SNPannot file, which contains SNP.pos, SNP.ref,SNP.alt; The rowAnnotations
      snpIndex += 1
      outFileSNPanno.writeLine(join([$rec.POS, $rec.REF[0],$rec_alt], sep="\t") )
      discard addViNode(barcodeTable = barcodeTable,  alleleCountTable = alleleCountTable,  scSpermSeq = scSpermSeq,
               outFileTotalCountMtx = outFileTotalCountMtx, outFileAltCountMtx  = outFileAltCountMtx, nnsize = nnsize,
               mindp = mindp, maxdp = maxdp, thetaRef = thetaRef, thetaAlt = thetaAlt, snpIndex = snpIndex, rec=rec,
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
  let version = "0.3.2"
  var doc = format("""
  $version

  Usage:
      sgcocaller fmf [options] <BAM> <VCF> <barcodeFile> <out_prefix>
      sgcocaller xo [options] <BAM> <VCF> <barcodeFile> <out_prefix>

Arguments:

  <BAM> the read alignment file with records of single-cell DNA reads
  
  <VCF> the variant call file with records of SNPs

  <barcodeFile> the text file containing the list of cell barcodes

  <out_prefix>  the prefix of output files


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
  -h --help  show help


  Examples
      ./sgcocaller fmf --threads 4 AAAGTAGCACGTCTCT-1.raw.bam AAAGTAGCACGTCTCT-1.raw.bam.dp3.alt.vcf.gz barcodeFile.tsv ./percell/cellfragment.fmf
      ./sgcocaller xo --threads 4 AAAGTAGCACGTCTCT-1.raw.bam AAAGTAGCACGTCTCT-1.raw.bam.dp3.alt.vcf.gz barcodeFile.tsv ./percell/ccsnp
""" % ["version", version])

  let args = docopt(doc, version=version)
  
  var
    threads:int
    vcff:string
    barcodeFile:string
    bamfile:string
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
  echo $args
  threads = parse_int($args["--threads"])
  barcodeTag = $args["--barcodeTag"]
  mindp = parse_int($args["--minDP"])
  maxdp = parse_int($args["--maxDP"])
  maxtotal = parse_int($args["--maxTotalDP"])
  mintotal = parse_int($args["--minTotalDP"])
  minsnpdepth = parse_int($args["--minSNPdepth"])

  bamfile = $args["<BAM>"]
  barcodeFile = $args["<barcodeFile>"]
  out_dir = $args["<out_prefix>"]
  vcff = $args["<VCF>"]
  mapq = parse_int($args["--minMAPQ"])
  minbsq = parse_int($args["--baseq"])
  thetaRef = parse_float($args["--thetaREF"])
  thetaAlt = parse_float($args["--thetaALT"])
  cmPmb = parse_float($args["--cmPmb"])
  
  var 
    ibam:Bam
    ivcf:VCF
    s_Chrs: seq[string]
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
  var barcodeTable =  newOrderedTable[string,int](initialSize = 1024)
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  ## initiate the table with CB as keys, allele counts (ref object) as elements
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if $kstr.s[0] == "#":
      continue
    var v = $kstr.s
    discard barcodeTable.hasKeyOrPut(v, 0)
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
                     maxtotal, mindp, maxdp, thetaREF, thetaALT, cmPmb,s_Chrs,barcodeTag)

  ibam.close()
  ivcf.close()