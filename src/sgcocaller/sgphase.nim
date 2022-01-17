# phase from genotype matrices generated from sgcocaller gtMtx

import findTemplateCell
import inferMissingSnps
import utils
import streams
import strutils
import sequtils

proc writePhasedSnpAnnot*(fullGeno:seq[BinaryGeno], snpAnnotFileStream:FileStream, phasedSnpAnnotFileStream:FileStream): int = 
  var snpindex = 0
  var header: string
  var snp_rec: string
  header = snpAnnotFileStream.readLine()
  phasedSnpAnnotFileStream.writeLine(join([header,"Phase"],sep = "\t"))
  while not snpAnnotFileStream.atEnd():
    snp_rec = snpAnnotFileStream.readLine()
    phasedSnpAnnotFileStream.writeLine(join([snp_rec,$fullGeno[snpindex]], sep  = "\t") )
    snpindex.inc
  if snpindex != (fullGeno.len):
    quit "snpAnnot.txt does not have the same number of rows with gtMtx"
  return 0 
 
proc sgphase*(mtxFile: string, snpAnnotFile:string, phasedSnpAnnotFile:string, diagnosticDataframeFile:string, templateCell:int, maxExpand = 1000, posteriorProbMin = 0.99,maxDissim:float): int = 
  var nsnps,ncells,totalEntries:int
  var gtMtx:seq[seq[BinaryGeno]]  
  var currentEntry:seq[int]
  var currentEntrySeq:seq[string]
  var gtMtxFileStream,ddframFileStream,snpAnnotFileStream,phasedSnpAnnotFileStream:FileStream
  #let phasedSnpAnnotFile = phaseOutdir & "phased_snpAnnot.txt"
  #let diagnosticDataframeFile = phaseOutdir & "cellGenoVersusTemplate.txt"
  ## i, j are 1-based, entries with value 1 or 2
  try: 
    gtMtxFileStream = openFileStream(mtxFile, fmRead)
    ddframFileStream = openFileStream(diagnosticDataframeFile, fmWrite)
    snpAnnotFileStream = openFileStream(snpAnnotFile, fmRead)
    phasedSnpAnnotFileStream = openFileStream(phasedSnpAnnotFile, fmWrite)
  except:
    stderr.write getCurrentExceptionMsg()
  # %%MatrixMarket matrix coordinate integer general
  discard gtMtxFileStream.readLine()
  #N, i,j
  currentEntrySeq = gtMtxFileStream.readLine().splitWhitespace()
  currentEntry = map(currentEntrySeq, proc(x: string): int = parseInt(x))
  nsnps = currentEntry[0]
  echo "nsnps " & $nsnps
  ncells = currentEntry[1]
  echo "ncells " & $ncells
  totalEntries = currentEntry[2]
  echo "totalEntries " & $totalEntries
  ## gtMtx is cell by Snp format
  gtMtx = newSeqWith(ncells,newSeq[BinaryGeno](nsnps))
  discard readGtMtxToSeq(gtMtxFileStream,gtMtx)
  var tempcell: int
  if templateCell != -1:
    tempcell = templateCell
  else:
    tempcell = selectTemplateCell(gtMtx = gtMtx, nPairs =3, maxDissim = maxDissim)
  echo "template cell is " & $tempcell
  var tempcellGeno = gtMtx[tempcell]
  var fullGeno = inferSnpGeno(tempcellGeno, gtMtx,posterior_thresh = posteriorProbMin, runCorrection = false, maxExpand = maxExpand)
  echo "inferSnoGeno done"
  # echo "second pass"
  # fullGeno = inferSnpGeno(fullGeno, gtMtx)
  # echo "second pass done"
  # generate the txt for diagnositc plot
  echo "run correction"
  fullGeno = inferSnpGeno(fullGeno, gtMtx,posterior_thresh = posteriorProbMin, runCorrection = true, maxExpand = maxExpand)
  echo "run correction done"
  var selectedCells:seq[int]
  if ncells < 10:
    selectedCells = (0..(ncells-1)).toSeq
  else:
    selectedCells = (0..9).toSeq
  var headliner = "templateGeno"
  var writeOut = ""
  for c in selectedCells:
    headliner = headliner & " cell" & $c
  ddframFileStream.writeLine(headliner)
  let subsetGtMtx = map(selectedCells, proc(x: int): seq[BinaryGeno] = gtMtx[x] )
  for k,g in fullGeno:
    if g == gUnknown: continue
    else:
      writeOut = $(parseInt($g) + 1) 
      for cellg in sliceColumn(subsetGtMtx, k):
        writeOut = writeOut & " " & $(parseInt($cellg)+1)
    ddframFileStream.writeLine(writeOut)
  discard writePhasedSnpAnnot(fullGeno, snpAnnotFileStream, phasedSnpAnnotFileStream)
  for fs in [gtMtxFileStream,ddframFileStream,snpAnnotFileStream,phasedSnpAnnotFileStream]:
    fs.close()
