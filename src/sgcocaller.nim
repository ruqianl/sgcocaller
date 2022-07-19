
import os
import docopt
import strutils
import hts
import tables
import sequtils
import sgcocaller/utils
import math
import streams
import sgcocaller/graph
import sgcocaller/findPath
import sgcocaller/getGtMtx
import sgcocaller/sgphase
import sgcocaller/writeVCF
import sgcocaller/correctPhase
import sgcocaller/sgcocaller_sxo
import times

let initProb:array[stateRef..stateAlt, float]=[0.5,0.5]

proc sgcocaller(threads:int, ivcf:VCF, barcodeTable:TableRef,
                ibam:Bam, out_dir:string, mapq:int, 
                minbsq:int, mintotal:int, maxtotal:int, mindp:int, maxdp:int, 
                thetaREF:float, thetaALT:float, cmPmb:float,s_Chrs:seq,barcodeTag:string,phased:bool): int =
  
  var bulkBam = false
  if barcodeTable.len == 1 and barcodeTable.hasKey("bulk"):
    echo $now() & " running in bulk mode and all reads are regarded as from one sample/cell"
    bulkBam = true
  else:
    echo $now() & " running sgcoaller xo for " & $barcodeTable.len & " single cells"

  ## iterate through each selected chromosome
  for chrom in s_Chrs:
    var snpIndex, nnsize = 0
    ## number of non zeros
    var scSpermSeq = newTable[int,SpermViNodes]()
    ## matches with the order in barcodeTable
    ## scSpermSeq.setLen(barcodeTable.len)
    var outFileSNPanno,outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx,viSegmentInfo:FileStream
    let sparseMatrixHeader = "%%MatrixMarket matrix coordinate integer general"
    try:
      outFileSNPanno = openFileStream(out_dir & chrom & "_snpAnnot.txt", fmWrite)
      outFileTotalCountMtx = openFileStream(out_dir & chrom & "_totalCount.mtx", fmWrite)
      outFileAltCountMtx = openFileStream(out_dir & chrom & "_altCount.mtx", fmWrite)
      outFileVStateMtx = openFileStream(out_dir & chrom & "_vi.mtx", fmWrite)
      viSegmentInfo = openFileStream(out_dir & chrom & "_viSegInfo.txt", fmWrite)
    except:
      stderr.write getCurrentExceptionMsg()    
    ## write headers to those mtx files and the first line place holder for total_row total_column total_entry    
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
                                        maxTotalReads = maxtotal, minTotalReads = mintotal,bulkBam = bulkBam, barcodeTag = barcodeTag,
                                        startPos = rec.POS.cint-1, stopPos=rec.POS.cint,rec_alt =  rec_alt, rec_ref = rec_ref)
      if alleleCountTable.len==0: continue

      ## add to snpAnnoSeq, later write to SNPannot file, which contains SNP.pos, SNP.ref,SNP.alt; The rowAnnotations
      snpIndex += 1
      outFileSNPanno.writeLine(join([$rec.POS, $rec_ref, $rec_alt], sep="\t") )
      discard addViNode(barcodeTable = barcodeTable,  alleleCountTable = alleleCountTable,  scSpermSeq = scSpermSeq,
               outFileTotalCountMtx = outFileTotalCountMtx, outFileAltCountMtx  = outFileAltCountMtx, nnsize = nnsize,
               mindp = mindp, maxdp = maxdp, thetaRef = thetaRef, thetaAlt = thetaAlt, snpIndex = snpIndex, rec_pos=int(rec.POS),
               initProb = initProb, cmPmb = cmPmb)
    discard pathTrackBack(scSpermSeq = scSpermSeq, thetaRef = thetaRef,thetaAlt=thetaAlt,cmPmb = cmPmb,outFileVStateMtx = outFileVStateMtx,
                          viSegmentInfo = viSegmentInfo )

    for outFileStream in [outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx]:
      outFileStream.setPosition(49)
      outFileStream.write($snpIndex & " " & $barcodeTable.len & " " & $nnsize) 
    for outFileStream in [outFileTotalCountMtx,outFileAltCountMtx,outFileVStateMtx,outFileSNPanno,viSegmentInfo]:
      outFileStream.close()
  return 0

when(isMainModule):
  let version = "0.3.9"
  var doc = format("""
  $version
  Usage:
      sgcocaller autophase [options] <BAM> <VCF> <barcodeFile> <out_prefix>
      sgcocaller phase [options] <BAM> <VCF> <barcodeFile> <out_prefix> 
      sgcocaller swphase [options] <gtMtxFile> <phasedSnpAnnotFile> <referenceVCF> <out_prefix> 
      sgcocaller sxo [options] <SNPPhaseFile> <phaseOutputPrefix> <barcodeFile> <out_prefix>
      sgcocaller xo [options] <BAM> <VCF> <barcodeFile> <out_prefix>
      

Arguments:

  <BAM> the read alignment file with records of single-cell DNA reads
  
  <VCF> the variant call file with records of SNPs

  <barcodeFile> the text file containing the list of cell barcodes

  <out_prefix>  the prefix of output files

  <out_vcf> the output vcf aftering phasing blocks in hapfile
  
  <gtMtxFile> the output gtMtx.mtx file from running sgcocaller phase

  <phasedSnpAnnotFile>  the output phased_snpAnnot.txt from running sgcocaller phase

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
  --minSNPdepth <minSNPdepth>  the minimum depth of cell coverage for a SNP to be includes in generated genotype matrix file [default: 1]
  --thetaREF <thetaREF>  the theta for the binomial distribution conditioning on hidden state being REF [default: 0.1]
  --thetaALT <thetaALT>  the theta for the binomial distribution conditioning on hidden state being ALT [default: 0.9]
  --cmPmb <cmPmb>  the average centiMorgan distances per megabases default 0.1 cm per Mb [default: 0.1]
  --phased  the input VCF for calling crossovers contains the phased GT of heterozygous SNPs
  --outvcf  generate the output in vcf format (sgcocaller phase)  
  --templateCell <templateCell>  the cell's genotype to be used a template cell, as the cell's index (0-starting) in the barcode file, default as not supplied [default: -1]
  --maxDissim <maxDissim>  the maximum dissimilarity for a pair of cell to be selected as potential template cells due to not having crossovers in either cell [default: 0.0099]
  --maxExpand <maxExpand>  the maximum number of iterations to look for locally coexisting positions for inferring missing SNPs in template haplotype sequence [default: 1000]
  --posteriorProbMin <posteriorProbMin>  the min posterior probability for inferring missing SNPs [default: 0.99]
  --lookBeyondSnps <lookBeyondSnps>  the number of local SNPs to use when finding switch positions [default: 25]
  --minSwitchScore <minSwitchScore>  the minimum switch score for a site to be identified as having a switch error in the inferred haplotype  [default: 50.0]
  --minPositiveSwitchScores <minPositiveSwitchScores>  the min number of continuing SNPs with positive switch scores to do switch error correction [default: 8]  
  --binSize <binSize>  the size of SNP bins for scanning swith errors, users are recommended to increase this option when SNP density is high. [default: 2000]
  --stepSize <stepSize>  the move step size used in combination with --binSize. [default: 200]
  --dissimThresh <dissimThresh>  the threshold used on the allele concordance ratio for determining if a SNP bin contains a crossover. [default: 0.0099]
  --batchSize <batchSize>  the number of cells to process in one batch when running sxo. This option is only needed when the memory is limited. 
  --notSortMtx  do not sort the output mtx. 
  --maxUseNcells <maxUseNcells>  the number of cells to use for calculating switch scores. All cells are used if not set
  -h --help  show help


  Examples
      ./sgcocaller autophase possorted_bam.bam hetSNPs.vcf.gz barcodeFile.tsv phaseOutputPrefix
      ./sgcocaller phase possorted_bam.bam hetSNPs.vcf.gz barcodeFile.tsv phaseOutputPrefix
      ./sgcocaller xo --threads 4 possorted_bam.bam phased_hetSNPs.vcf.gz barcodeFile.tsv ./percell/ccsnp
      ./sgcocaller sxo snp_phase.txt phaseOutputPrefix barcodeFile.tsv ./percell/ccsnp

""" % ["version", version])

  let args = docopt(doc, version=version)
  
  var
    threads,mapq,minbsq,mintotal,maxtotal,mindp,maxdp,minsnpdepth,maxExpand,lookBeyondSnps,minPositiveSwitchScores,binSize,stepSize,batchSize,maxUseNcells:int
    thetaREF,thetaALT,cmPmb,posteriorProbMin,maxDissim,minSwitchScore,dissimThresh:float
    barcodeTag="CB"
    out_dir,selectedChrs,barcodeFile,bamfile,vcff:string
    vcfGtPhased,notSortMtx = false
    s_Chrs: seq[string]
    outFileTotalCountMtxFile, outFileAltCountMtxFile: string

  echo $now() & $args
  echo $now() & " running sgcocaller " & $version
  threads = parse_int($args["--threads"])
  barcodeTag = $args["--barcodeTag"]
  mindp = parse_int($args["--minDP"])
  maxdp = parse_int($args["--maxDP"])
  maxtotal = parse_int($args["--maxTotalDP"])
  mintotal = parse_int($args["--minTotalDP"])
  minsnpdepth = parse_int($args["--minSNPdepth"])
  mapq = parse_int($args["--minMAPQ"])
  minbsq = parse_int($args["--baseq"])
  maxExpand = parse_int($args["--maxExpand"])
  lookBeyondSnps = parse_int($args["--lookBeyondSnps"])
  binSize = parse_int($args["--binSize"])
  stepSize = parse_int($args["--stepSize"])
  dissimThresh = parse_float($args["--dissimThresh"])
  minPositiveSwitchScores =  parse_int($args["--minPositiveSwitchScores"])
  thetaRef = parse_float($args["--thetaREF"])
  thetaAlt = parse_float($args["--thetaALT"])
  posteriorProbMin = parse_float($args["--posteriorProbMin"])
  maxDissim = parse_float($args["--maxDissim"])
  minSwitchScore = parse_float($args["--minSwitchScore"])
  if $args["--maxUseNcells"] != "nil": maxUseNcells = parse_int($args["--maxUseNcells"])
  cmPmb = parse_float($args["--cmPmb"])
  vcfGtPhased = parse_bool($args["--phased"])
  notSortMtx = parse_bool($args["--notSortMtx"])
  out_dir = $args["<out_prefix>"]
  if (out_dir[^1] != '_') and (out_dir[^1] != '/') : out_dir = out_dir & "_"
  if $args["--batchSize"] != "nil":
    batchSize = parse_int($args["--batchSize"])
  else:
    batchSize = 0
  var run_phase, run_swphase, run_xo, run_sxo, run_auto_phase = false
  run_phase = args["phase"]
  run_swphase = args["swphase"]
  run_xo = args["xo"]
  run_sxo = args["sxo"]
  run_auto_phase = args["autophase"]
  if run_auto_phase:
    run_phase = true
    run_swphase = true
  if run_phase or run_xo or run_sxo:
    barcodeFile = $args["<barcodeFile>"]
    var hf = hts.hts_open(cstring(barcodeFile), "r")
    var barcodeTable =  newTable[string,int](initialSize = 1024)
    var kstr: hts.kstring_t
    kstr.l = 0
    kstr.m = 0
    kstr.s = nil
    var ithSperm = 0
    ## initiate the table with CB as keys, gamete index as elements
    while hts_getline(hf, cint(10), addr kstr) > 0:
      if $kstr.s[0] == "#": continue
      var v = $kstr.s
      discard barcodeTable.hasKeyOrPut(v, ithSperm)
      ithSperm.inc
    discard hf.hts_close()
    if run_phase or run_xo:
      var 
        ibam:Bam
        ivcf:VCF
      bamfile = $args["<BAM>"]
      vcff = $args["<VCF>"]
      if($args["--chrom"] != "nil"):
        selectedChrs = $args["--chrom"]
        s_Chrs = selectedChrs.split(',')
      else:
        ## find all chroms from headers of vcf file (Contigs)
        let contigs = ivcf.contigs 
        s_Chrs = map(contigs, proc(x: Contig): string = x.name)      
      if not open(ibam, bamfile, threads=threads, index = true):
          quit "couldn't open input bam"
      if not open(ivcf, vcff, threads=threads):
          quit "couldn't open: vcf file"
      if run_phase:
        echo $now() & " running phase"
        var templateCell = parse_int($args["--templateCell"])     
        var outGtMtxFile, outTotalCountMtxFile,outAltCountMtxFile, outSnpAnnotFile, outphasedSnpAnnotFile,outdiagnosticDataframeFile : string
        var outGtMtx, outSnpAnnot, outTotalCountMtx, outAltCountMtx:FileStream
        for chrom in s_Chrs:
          echo $now() & " generating genotype sparse matrix file to " & out_dir & "for chr: " & chrom
          outGtMtxFile = out_dir & chrom & "_gtMtx.mtx"
          outTotalCountMtxFile = out_dir & chrom &  "_totalMtx.mtx"
          outAltCountMtxFile = out_dir & chrom &  "_altMtx.mtx"
          outSnpAnnotFile = out_dir & chrom &  "_snpAnnot.txt"
          outphasedSnpAnnotFile = out_dir & chrom &  "_phased_snpAnnot.txt"
          outdiagnosticDataframeFile = out_dir & chrom &  "_cellGenoVersusTemplate.txt"
          try: 
            outGtMtx = openFileStream(outGtMtxFile, fmWrite)
            outSnpAnnot = openFileStream(outSnpAnnotFile, fmWrite)
            outTotalCountMtx = openFileStream(outTotalCountMtxFile, fmWrite)
            outAltCountMtx = openFileStream(outAltCountMtxFile, fmWrite)
          except:
            stderr.write getCurrentExceptionMsg()
          discard getGtMtx(ibam = ibam, ivcf = ivcf, barcodeTable = barcodeTable, outGtMtx = outGtMtx, outTotalCountMtx = outTotalCountMtx ,outAltCountMtx = outAltCountMtx, outSnpAnnot = outSnpAnnot, maxTotalReads = maxtotal, minTotalReads = mintotal, mapq = mapq, minbsq = minbsq, minCellDp = mindp,maxCellDp = maxdp,barcodeTag = barcodeTag, minsnpdepth = minsnpdepth, chrom = chrom)                    
          if not notSortMtx:
            for mtxFile in [outGtMtxFile, outTotalCountMtxFile, outAltCountMtxFile]:
              var imtx = readMtx(mtx_file = mtxFile)
              discard sortWriteMtx(imtx, mtx_file = mtxFile) 
          echo $now() & "running sgcocaller phase from single cell genotype matrix (of one chromosome) " & outGtMtxFile
          echo $now() & "using the " & outSnpAnnotFile & " for generating phased haplotypes"
          discard sgphase(mtxFile = outGtMtxFile, snpAnnotFile = outSnpAnnotFile, phasedSnpAnnotFile = outphasedSnpAnnotFile,
                          diagnosticDataframeFile = outdiagnosticDataframeFile, templateCell = templateCell, 
                          maxExpand = maxExpand, posteriorProbMin = posteriorProbMin,maxDissim = maxDissim)
          if parse_bool($args["--outvcf"]) and (not run_swphase):
            var outvcfFile = out_dir & chrom &  "_phased_snpAnnot.vcf.gz"  
            echo $now() & " write phased haplotypes to VCF " & outvcfFile
            discard writePhaseToVCF(vcff, outvcfFile, outphasedSnpAnnotFile,threads = threads)    
          if run_swphase:
            echo $now() & " running swphase to find switch spots and generate corrected phase"
            let switchedPhasedAnnotFile = out_dir & chrom & "_corrected_phased_snpAnnot.txt"
            let switchScoreFile = out_dir & chrom & "_switch_score.txt" 
            let switchedPhasedAnnotVcfFile  = out_dir & chrom & "_corrected_phased_snpAnnot.vcf.gz"
            echo $now() & "scanning with swphase to find switch spots and generate corrected phase"
            discard correctPhase(outGtMtxFile, outphasedSnpAnnotFile, switchedPhasedAnnotFile, switchScoreFile, lookBeyondSnps = lookBeyondSnps, minSwitchScore = minSwitchScore,
                                 minPositiveSwitchScores = minPositiveSwitchScores, binSize = binSize, stepSize = stepSize,dissimThresh = dissimThresh, maxUseNcells = maxUseNcells )
            echo $now() & " write corrected phased haplotypes to " & switchedPhasedAnnotVcfFile
            discard writePhaseToVCF(vcff, switchedPhasedAnnotVcfFile, switchedPhasedAnnotFile, add_header_string = """##sgcocaller_v0.3.9=swphase""",threads = threads, chrom = chrom)    
      elif run_xo:
        echo $now() & " running crossover calling from VCF and BAM file\n"
        echo $now() & " running for these chromosomes: " & $s_Chrs
        discard sgcocaller(threads, ivcf, barcodeTable, ibam, out_dir, mapq, minbsq, mintotal, maxtotal, mindp, maxdp, thetaREF, thetaALT, cmPmb, s_Chrs, barcodeTag, phased = vcfGtPhased)
        ## sort entries in mtx files 
        if not notSortMtx:
          for chrom in s_Chrs:
            outFileTotalCountMtxFile = out_dir & chrom & "_totalCount.mtx"
            outFileAltCountMtxFile = out_dir & chrom & "_altCount.mtx"
            for mtxFile in [outFileTotalCountMtxFile, outFileAltCountMtxFile]:
              var imtx = readMtx(mtx_file = mtxFile)
              discard sortWriteMtx(imtx, mtx_file = mtxFile)
      ibam.close()
      ivcf.close()
    elif run_sxo:
      if($args["--chrom"] != "nil"):
        selectedChrs = $args["--chrom"]
        s_Chrs = selectedChrs.split(',')
      else:
        quit "chromosomes need to be supplied via --chrom. Supply multiple chromosomes (will be processed sequencially) using comma as separator"
      echo $now() & " running crossover calling from genotype and allele count matrices generated from sgcocaller phase/swphase "
      echo $now() & " running for these chromosomes : " & $s_Chrs
      let phase_dir = $args["<phaseOutputPrefix>"]
      let phasedSnpAnnotFile = $args["<SNPPhaseFile>"]
      if (batchSize == 0) or (batchSize >= barcodeTable.len):
        batchSize = barcodeTable.len        
      discard sgcocallerSXO(barcodeTable = barcodeTable, phase_dir = phase_dir, out_dir = out_dir,
                            thetaREF = thetaREF, thetaALT = thetaALT, cmPmb = cmPmb,s_Chrs =s_Chrs,initProb = initProb,
                            phasedSnpAnnotFileName = phasedSnpAnnotFile,batchSize=batchSize)
      if not notSortMtx:
        for chrom in s_Chrs:
          outFileTotalCountMtxFile = out_dir & chrom & "_totalCount.mtx"
          outFileAltCountMtxFile = out_dir & chrom & "_altCount.mtx"
          for mtxFile in [outFileTotalCountMtxFile, outFileAltCountMtxFile]:
            var imtx = readMtx(mtx_file = mtxFile)
            discard sortWriteMtx(imtx, mtx_file = mtxFile)
  elif run_swphase:
    echo $now() & " running swphase to find switch spots and generate corrected phase"
    let gtMtxFile =  $args["<gtMtxFile>"]
    let phasedSnpAnnotFile = $args["<phasedSnpAnnotFile>"]
    let switchedPhasedAnnotFile = out_dir & "corrected_phased_snpAnnot.txt"
    let switchScoreFile = out_dir & "switch_score.txt" 
    let switchedPhasedAnnotVcfFile  = out_dir & "corrected_phased_snpAnnot.vcf.gz"
    echo $now() & " scanning with swphase to find switch spots and generate corrected phase"
    discard correctPhase(gtMtxFile,phasedSnpAnnotFile,switchedPhasedAnnotFile,switchScoreFile,lookBeyondSnps = lookBeyondSnps,minSwitchScore = minSwitchScore, 
                         minPositiveSwitchScores = minPositiveSwitchScores, binSize = binSize, stepSize = stepSize, dissimThresh =dissimThresh,maxUseNcells = maxUseNcells)
    if ($args["--chrom"] == "nil"):
      echo $now() & " assuming supplied VCF only contains the SNPs for the relevant chromosome. If this is not the case, use the --chrom option"
    echo $now() & " write corrected phased haplotypes to " & switchedPhasedAnnotVcfFile
    discard writePhaseToVCF($args["<referenceVCF>"], switchedPhasedAnnotVcfFile, switchedPhasedAnnotFile, add_header_string = """##sgcocaller_v0.3.9=swphase""",threads = threads, chrom = $args["--chrom"])




 
