from distributions/rmath import dbinom
import math
import tables
import hts

type 
  allele_expr* = ref object
    cref*: int
    calt*: int
  ViState* = enum
    stateRef, stateAlt, stateN
  ViNode* = object
    pos*: int
    cRef*: int
    cAlt*: int
    pathScore*: array[stateRef..stateAlt, float]
    pathState*: array[stateRef..stateAlt, ViState]      
    state*: ViState  
  GtNode* = ref object
    chrom*: string
    variantIndex*: int
    alleles*: string
    genotype*: int
  Fragment* = ref object
    cellbarcode*: string
    counter*: int
    node1*: GtNode
    node2*: GtNode
proc get_base_offset*(position:int, align: Record): int =
  var 
    off = align.start.int
    qoff = 0
    roff_only = 0
    base_offset = 0
    over = 0
  for event in align.cigar:
    var cons = event.consumes
    if cons.query:
      qoff+=event.len ## the offs to the query sequences
    if cons.reference:
      off += event.len ## the offs to the reference sequences
      if not cons.query:
        roff_only += event.len
    if off <= position:
      continue
    over = off - position
    # get the base 
    base_offset = qoff - over-1
    break 
  return base_offset   
proc add_allele*(g:GtNode, alleleBinGeno:int):string = g.alleles & $alleleBinGeno

proc inc_count_cref*(a:allele_expr) = inc(a.cref)
proc inc_count_calt*(a:allele_expr) = inc(a.calt)

proc getEmission*(thetaRef=0.1,thetaAlt=0.9,
                  cRef:int,cAlt:int): array[stateRef..stateAlt, float] =
  
  var emissionScore = [dbinom(x=float(cAlt),size=(cRef+cAlt),prob=thetaRef,log=true),
                        dbinom(x=float(cAlt),size=(cRef+cAlt),prob=thetaAlt,log=true)]
  return emissionScore

proc getTrans*(pos1:int64,pos2:int64,cmPmb=0.1): float = 
  if pos2<pos1:
    quit "Wrong order of snps"
  var rec = 1-math.exp(-float(pos2-pos1)*1e-8*cmPmb) # for autosomes 1*cmPbm cM per Mb 
  return rec
                   #   
proc countAllele*(ibam:Bam, maxTotalReads:int,
                  minTotalReads:int,
                  chrom:string, mapq: int,
                  barcodeTable:TableRef,minbsq:int,bulkBam:bool,barcodeTag:string,startPos:int, stopPos:int,
                  rec_alt:char, rec_ref:char): Table[string,allele_expr] =
  var alleleCountTable = initTable[string,allele_expr]()
  #var rec_alt:char
  var total_reads = 0
  var base_off: int
  var base: char
  #rec_alt = rec.ALT[0][0]
  for aln in ibam.query(chrom = chrom,start = startPos, stop = stopPos):
    var cbt = tag[string](aln, barcodeTag)
    var currentCB: string
    if cbt.isNone:
      if not bulkBam:
        continue  
      else:
        currentCB = "bulk"
    else:
      currentCB = cbt.get  
    if aln.flag.unmapped or aln.mapping_quality.cint < mapq or aln.flag.dup: continue
    ## skip unmapped, duplicated, mapq low reads or aln.flag.secondary or aln.flag.supplementary
    ## if not aln.flag.proper_pair: continue
    if not barcodeTable.hasKey(currentCB): continue
    base_off = get_base_offset(position = stopPos, align = aln) 
    base = aln.base_at(base_off)
    if aln.base_quality_at(base_off).cint < minbsq: 
      continue
    total_reads+=1
    if alleleCountTable.hasKey(currentCB):
      if aln.base_at(base_off) == rec_ref: 
        alleleCountTable[currentCB].inc_count_cref
        continue
      if aln.base_at(base_off) == rec_alt: 
        alleleCountTable[currentCB].inc_count_calt
        continue
    else:
      var new_snp = allele_expr(cref:0,calt:0)
      alleleCountTable[currentCB] = new_snp
      if aln.base_at(base_off) == rec_ref: 
        alleleCountTable[currentCB].inc_count_cref
        continue
      if aln.base_at(base_off) == rec_alt: 
        alleleCountTable[currentCB].inc_count_calt
        continue
  if total_reads > maxTotalReads or total_reads < minTotalReads:
    return initTable[string,allele_expr]()
  return alleleCountTable