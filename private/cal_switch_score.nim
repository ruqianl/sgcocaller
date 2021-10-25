## calcualte switch score for snps positions that are high risk
import utils
import sequtils
import math
import streams
# bins of SNP indexes
import strutils

let binSize = 1000
let movingStep = 500
let dissimThresh = 0.0099 
let lookBeyondSnps = 25
let debug = false
let switchPrior = 0.5
type 
  switchScoreTuple = tuple
    switch_scores: seq[float]
    switch_scores_snpIndex: seq[int]
    
# reduce the search space
# return seq of snp indexes
# low risk bins are bins that do not have crossover happening in majority of the cells. 

## rough crossover estimates: return binary results, 0 means no crossovers, 1 means there is.
## simply comparing the dissimilary btw template geno seq with cell's geno seq for the SNPs in the bin

## cell_geno snp by cell
proc hasCrossover(templ_geno:seq[BinaryGeno], cell_geno: seq[seq[BinaryGeno]]): bool = 
  let ncells = cell_geno[0].len
  var ncompared = newSeq[int](ncells)
  var nmatch = newSeq[int](ncells)
  for snpi,g in templ_geno:
    for cellj in 0..(ncells-1):
      if g != gUnknown and cell_geno[snpi][cellj] != gUnknown:
        ncompared[cellj].inc
        if g == cell_geno[snpi][cellj]: nmatch[cellj].inc
  if debug: echo "ncompared " & $ncompared      
  let dissim = map(toSeq(0..(ncells-1)), proc(x:int):float =  nmatch[x]/ncompared[x])
  if debug: echo "dissim " & $dissim
  var ncrossover = map(dissim, proc(x:float): int = (int)(x > dissimThresh and x < (1-dissimThresh)))
  if debug: echo "ncrossover " & $ncrossover
  let nxcells = foldl(ncrossover, a + b)
  if nxcells < int(floor(0.5 * float(ncells))):
    if debug: echo "bin had no crossovers : " & $nxcells
    return false
  else:
    if debug: echo "bin had many crossovers : " & $nxcells
    return true

proc findHighRiskSnps(fullGeno:seq[BinaryGeno], gtMtxByCell:seq[seq[BinaryGeno]]): seq[int] = 
  #let ncells = gtMtxByCell[0].len
  let nSnps = gtMtxByCell.len
  let nBins = (int)(ceil(nSnps / (binSize - movingStep)))
  if debug: echo "nBins: " & $nBins
  let binIds = toSeq(0..(nBins-1))
  var snpPosStart,snpPosEnd: int
  var highRiskSnps: seq[int]
  for binI in binIds:
    if debug: echo "binI " & $binI
    snpPosStart = (binI)*(binSize - movingStep)
    if (snpPosStart + binSize) > (nSnps-1): 
      snpPosEnd = (nSnps-1)
    else:
      snpPosEnd = snpPosStart + binSize
    if (hasCrossover(templ_geno =  fullGeno[snpPosStart..snpPosEnd], 
                     cell_geno =  gtMtxByCell[snpPosStart..snpPosEnd] )):
      highRiskSnps = concat(highRiskSnps,toSeq(snpPosStart..snpPosEnd))

    
  # add a last bin as from end of chrom with length binSize
  let lastBinStart = nSnps - binSize
  let lastBinEnd = nSnps - 1
  if (hasCrossover(templ_geno =  fullGeno[snpPosStart..snpPosEnd], 
                     cell_geno =  gtMtxByCell[snpPosStart..snpPosEnd])):
    highRiskSnps = concat(highRiskSnps,toSeq(lastBinStart..lastBinEnd))

  highRiskSnps = deduplicate(highRiskSnps)
  return highRiskSnps
   
proc readPhasedSnpAnnot(phasedSnpAnnotFileStream: FileStream, nSnps:int): seq[BinaryGeno] =
  var fullGeno = newSeq[BinaryGeno](nSnps)
  var i = 0
  var currentEntrySeq: seq[string]
  while not phasedSnpAnnotFileStream.atEnd():
    currentEntrySeq = phasedSnpAnnotFileStream.readLine().splitWhitespace()
    fullGeno[i] = parseInt(currentEntrySeq[3]).toBinaryGeno(format = "01")
    i.inc
  if i != nSnps :
    quit "phased snpAnnot file does not have the same number of rows with gtMtx"
  return fullGeno

proc calculateProbslog10(hap:seq[BinaryGeno],cell_hap:seq[BinaryGeno], error_rate = 0.1) : float = 
  var nmatch,ncompare:int
  for i,g in hap:
    if (g != gUnknown) and (cell_hap[i] != gUnknown):
      ncompare.inc
      if (g == cell_hap[i]): nmatch.inc
  let nmis = ncompare - nmatch
  let prob = 0.5*(error_rate^nmis*(1-error_rate)^nmatch + error_rate^nmatch*(1-error_rate)^nmis)
  return(log10(prob))

proc countNmatch( hap:seq[BinaryGeno],cell_hap:seq[BinaryGeno],nmatch =true): int = 
  var nmatch,ncompare:int
  for i,g in hap:
    if (g != gUnknown) and (cell_hap[i] != gUnknown):
      ncompare.inc
      if (g == cell_hap[i]): nmatch.inc
  let nmis = ncompare - nmatch
  return nmatch

proc countNMis( hap:seq[BinaryGeno],cell_hap:seq[BinaryGeno],nmatch =true): int = 
  var nmatch,ncompare:int
  for i,g in hap:
    if (g != gUnknown) and (cell_hap[i] != gUnknown):
      ncompare.inc
      if (g == cell_hap[i]): nmatch.inc
  let nmis = ncompare - nmatch
  return nmis

proc binarySwitch(x: BinaryGeno):BinaryGeno = 
    if(x == gUnknown):  gUnknown 
    elif(x == gREF):gALT
    else: gREF

proc switchHap(hap: seq[BinaryGeno]): seq[BinaryGeno] = 
  map(hap, binarySwitch)

proc getIthCellHap(gtMtxByCell:seq[seq[BinaryGeno]],cellIndex:int):seq[BinaryGeno] = 
  map(gtMtxByCell, proc(y: seq[BinaryGeno]): BinaryGeno = y[cellIndex])

proc calSwitchScore(riskySnps: seq[int], gtMtxByCell:seq[seq[BinaryGeno]], fullGeno: seq[BinaryGeno]): switchScoreTuple = 
  var letfIndexStart,letfIndexEnd,rightIndexStart,rightIndexEnd,offset: int
  let nSnps = gtMtxByCell.len
  let ncells =  gtMtxByCell[0].len
  var hswitch,cell_hap,cell_full_hap,templ_geno,templ_geno_left,templ_geno_right: seq[BinaryGeno]
  var prob_switch, prob_nonsw: float
  var leftIndex,rightIndex, switch_score_snpIndex: seq[int]
  var swscore = newSeq[float](nSnps)
  for rsnpi in riskySnps:
    prob_switch = 0.0
    prob_nonsw = 0.0
    if fullGeno[rsnpi] == gUnknown: continue
    for celli in 0..(ncells-1):
      # if celli == 2: continue
      leftIndex = newSeq[int]()
      rightIndex = newSeq[int]()
      offset = 1
      cell_full_hap = getIthCellHap(gtMtxByCell, celli)
      if(cell_full_hap[rsnpi] != gUnknown and fullGeno[rsnpi] != gUnknown ):
        rightIndex.add(rsnpi)
      ## need to find the nearest coexisting N snps (N = lookBeyondSnps*2)
      offset = 1
      while(leftIndex.len < lookBeyondSnps and ((rsnpi-offset) > 0)):
        if((cell_full_hap[rsnpi-offset] != gUnknown) and (fullGeno[rsnpi-offset] != gUnknown)):
          leftIndex.add(rsnpi-offset)
        offset.inc
      offset = 1
      while(rightIndex.len < lookBeyondSnps and ((rsnpi+offset) < nSnps)):
        if((cell_full_hap[rsnpi+offset] != gUnknown) and (fullGeno[rsnpi+offset] != gUnknown)):
          rightIndex.add(rsnpi+offset)
        offset.inc
      cell_hap = concat(map(concat(leftIndex,rightIndex), proc(x:int):BinaryGeno = cell_full_hap[x])) 
      templ_geno_left = map(leftIndex, proc(x:int):BinaryGeno = fullGeno[x])
      templ_geno_right = map(rightIndex, proc(x:int):BinaryGeno = fullGeno[x])
      templ_geno = concat(templ_geno_left,templ_geno_right)
      hswitch = concat(templ_geno_left,switchHap(templ_geno_right))
      # if rsnpi == 16933:
      #   echo "cell " & $celli
      #   echo "cell_hap " & $cell_hap
      #   echo "templ_geno " & $templ_geno
      #   echo "hswitch " & $hswitch
      #   echo "prev prob nosw " & $prob_nonsw
      #   echo "prev prob_switch " & $prob_switch   
      #   echo "nmatch temp_geno " & $(countNmatch(hap = templ_geno, cell_hap = cell_hap))
      #   echo "nmatch hswitch " & $(countNmatch(hap = hswitch, cell_hap = cell_hap))
      #   echo "nmis temp_geno " & $(countNMis(hap = templ_geno, cell_hap = cell_hap))
      #   echo "nmis hswitch " & $(countNMis(hap = hswitch, cell_hap = cell_hap))    
      prob_nonsw = prob_nonsw + calculateProbslog10(hap = templ_geno, cell_hap = cell_hap)
      prob_switch = prob_switch + calculateProbslog10(hap = hswitch, cell_hap = cell_hap)

    #if(prob_switch > prob_nonsw): echo "found positive switch score"
    swscore[rsnpi] =  log10(switchPrior) + prob_switch - log10(1-switchPrior) - prob_nonsw
    switch_score_snpIndex.add(rsnpi)
  return (swscore,switch_score_snpIndex)

## return SNP indexs to switch
## switchScore, a dense list of switch scores
proc findSwitchSites(switchScoreT: switchScoreTuple): seq[int] = 
  var sitesToSwitch = newSeq[int]()
  var positiveScores: seq[float]
  var positiveScoresIndex: seq[int]
  
  var inPositiveBlock = false
  for site, score in switchScoreT.switch_scores:
    if switchScoreT.switch_scores_snpIndex.find(site)>=0:
      if score <= 0.0 :
        if inPositiveBlock:
          inPositiveBlock = false
          if positiveScoresIndex.len > int(floor(lookBeyondSnps/2)):
            sitesToSwitch.add(positiveScoresIndex[maxIndex(positiveScores)])
          continue
      elif not inPositiveBlock:
        positiveScores = @[score]
        positiveScoresIndex = @[site]
        inPositiveBlock = true
        continue
      else:
        positiveScores.add(score)
        positiveScoresIndex.add(site)
  return sitesToSwitch
let vcf_file = "data/swapped/FVB_NJ.mgp.v5.snps.dbSNP142.homo.alt.modified.swapped.GT.chr1.vcf.gz"
let bam_file = "data/WC_CNV_43/WC_CNV_43.bam"
let barcode_file = "data/WC_CNV_43/WC_CNV_43_filteredBC.tsv"
var
  gtMtxFileStream,phaseSnpAnnotFileStream,switchScoreFileStream:FileStream


let s_Chrs = map(toSeq(1..19), proc(x:int):string = "chr" & $x)

var gtMtxByCell: seq[seq[BinaryGeno]]
var currentEntrySeq: seq[string]
var currentEntry: seq[int]
var nSnps,ncells:int
## sort entries in mtx files 
for chrom in s_Chrs:    
  let gtMtxFile = "data/WC_CNV_43/sgcocaller/gtMtx/" & chrom & "_gtMtx.mtx"
  let phaseSnpAnnotFile = "data/WC_CNV_43/sgcocaller/phase/" & chrom & "_phased_snpAnnot.txt"
  let switchScoreFile = "data/WC_CNV_43/sgcocaller/phase/" & chrom & "_switch_score.txt"  
  try: 
    gtMtxFileStream = openFileStream(gtMtxFile, fmRead)
    phaseSnpAnnotFileStream = openFileStream(phaseSnpAnnotFile, fmRead)
    switchScoreFileStream = openFileStream(switchScoreFile, fmWrite)
    
  except:
    stderr.write getCurrentExceptionMsg() 
  echo "reading gtMtx to by cell"
  discard gtMtxFileStream.readLine()
  discard phaseSnpAnnotFileStream.readLine()
  #N, i,j
  currentEntrySeq = gtMtxFileStream.readLine().splitWhitespace()
  currentEntry = map(currentEntrySeq, proc(x: string): int = parseInt(x))
  nSnps = currentEntry[0]
  echo "nSnps " & $nSnps
  ncells = currentEntry[1]
  echo "ncells " & $ncells

  echo "totalEntries " & $ currentEntry[2]
  ## gtMtx is cell by Snp format
  gtMtxByCell = newSeqWith(nSnps,newSeq[BinaryGeno](ncells))
  discard readGtMtxToSeq(mtxFileStream = gtMtxFileStream, gtMtx = gtMtxByCell, by_cell = true)
  var fullGeno = readPhasedSnpAnnot(phasedSnpAnnotFileStream = phaseSnpAnnotFileStream, nSnps = nSnps )

  var riskySnps = findHighRiskSnps(fullGeno = fullGeno, gtMtxByCell = gtMtxByCell)
  echo "riskySnps.len " & $riskySnps.len
  var switchScoresT = calSwitchScore(riskySnps = riskySnps, gtMtxByCell =gtMtxByCell, fullGeno = fullGeno)
  switchScoreFileStream.writeLine("#" & $findSwitchSites(switchScoresT))
  echo  $findSwitchSites(switchScoresT)
  for snpi,score in switchScoresT.switch_scores:
    if switchScoresT.switch_scores_snpIndex.find(snpi)>=0:
      switchScoreFileStream.writeLine($score)
    else: 
      switchScoreFileStream.writeLine("NA")
  for fs in [gtMtxFileStream,phaseSnpAnnotFileStream,switchScoreFileStream]:
    fs.close()
 
