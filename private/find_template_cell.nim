# find template gamete that serves as a initial haplotype for this individual

import streams
import sequtils
import strutils
import math
import algorithm
import utils

type 
  CellPairs = object
    cell1*: int
    cell2*: int
    coexistSnps*: int
    dissimilarity*: float
  
let
  minCoexit = 1000
  maxDissim = 0.0099

var nsnps,ncells,totalEntries:int

proc toBinaryGeno(x: int) : BinaryGeno = 
  if x == 1: return gREF
  if x == 2: return gALT
  quit "mtx geno is not 1 or 2"

proc readGtMtx*(mtxFileStream:FileStream, gtMtx:var seq[seq[BinaryGeno]]):int = 
  var currentEntrySeq:seq[int]
  var currentLine:seq[string]

  while not mtxFileStream.atEnd():
    currentLine = mtxFileStream.readLine().splitWhitespace()
    ## i j 1-based from file
    currentEntrySeq = map(currentLine, proc(x: string): int = parseInt(x))
    gtMtx[(currentEntrySeq[1]-1)][(currentEntrySeq[0]-1)] = currentEntrySeq[2].toBinaryGeno
  return 0

proc countCoexistMatch(seq1: seq[BinaryGeno], seq2: seq[BinaryGeno]): seq[int] = 
  if not (seq1.len == seq2.len):
    quit "cells' genotypes under comparison do not have the same length"
  var ncoexist, nmatch, cell1Snps, cell2Snps:int
  for i, g in seq1:
    if g != gUnknown: cell1Snps += 1
    if (seq2[i] != gUnknown): cell2Snps += 1
    if (g != gUnknown and seq2[i]!= gUnknown):
      ncoexist += 1
      if(g == seq2[i]):
        nmatch += 1
  return [cell1Snps,cell2Snps,nmatch,ncoexist].toSeq

proc findCellPairs(gtMtx:seq[seq[BinaryGeno]],nPairs = 3, nSnpsPercell: var seq[int]): seq[CellPairs] = 
  let ncells = gtMtx.len
  var genoMatch = @[0,0,0,0]
  var returnCp = newSeq[CellPairs](nPairs)
  var k = 0
  for celli in 0..(ncells-2):
    if k == nPairs: break
    for cellj in (celli+1)..(ncells-1):
      genoMatch = countCoexistMatch(seq1 = gtMtx[celli], seq2 = gtMtx[cellj])
      nSnpsPercell[celli] = genoMatch[0]
      nSnpsPercell[cellj] = genoMatch[1]     
      if genoMatch[3] >= minCoexit:
        if (genoMatch[2]/genoMatch[3] < maxDissim) or (genoMatch[2]/genoMatch[3] > (1 - maxDissim)):
          returnCp[k] = CellPairs(cell1: celli, cell2: cellj, coexistSnps:genoMatch[3], dissimilarity:(genoMatch[2]/genoMatch[3]))
          k += 1
      if k == nPairs: break
  returnCp


proc findCellBynSNPs*(nSnpsPercell:seq[int],q = 0.75): int = 
  # return the selected cells' index
  # not too many SNPs or too few
  var chosed:int
  var seqNsnps = nSnpsPercell
  sort(seqNsnps, system.cmp[int])
  int(floor(float(seqNsnps.len) * q))

proc selectTemplateCell*(gtMtx:seq[seq[BinaryGeno]], nPairs = 3): int = 
  var nSnpsCells = newSeq[int](gtMtx.len)
  var cellPairs = findCellPairs(gtMtx = gtMtx, nPairs = nPairs, nSnpsPercell = nSnpsCells)
  echo "cellPairs"
  echo cellPairs
  var selectedCell, icp: int
  if cellPairs[0].coexistSnps == 0 :
    # no good pairs found
    selectedCell = findCellBynSNPs(nSnpsPercell = nSnpsCells)
  else:
    let coExistSnps = map(cellPairs, proc(x: CellPairs) : int = x.coexistSnps)
    icp = maxIndex(coExistSnps)
  if nSnpsCells[cellPairs[icp].cell1] > nSnpsCells[cellPairs[icp].cell2]:
    echo "selected cell: " & $cellPairs[icp].cell1 & " nSnps: " & $nSnpsCells[cellPairs[icp].cell1]
    return cellPairs[icp].cell1
  else:
    echo "selected cell: " & $cellPairs[icp].cell2 & " nSnps: " & $nSnpsCells[cellPairs[icp].cell2]
    return cellPairs[icp].cell2

var gtMtx:seq[seq[BinaryGeno]]  
var currentEntry:seq[int]
var currentEntrySeq:seq[string]
var gtMtxFileStream:FileStream
## i, j are 1-based
let mtxFile = "/mnt/mcfiles/rlyu/Projects/sgcocaller/test_data/WC_CNV_chr1_nonzero.mtx"
echo mtxFile

try: 
  gtMtxFileStream = openFileStream(mtx_file, fmRead)
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


gtMtx = newSeqWith(ncells,newSeq[BinaryGeno](nsnps))

discard readGtMtx(gtMtxFileStream,gtMtx)
#var colIndex = 4
#echo gtMtx[0][27..40]

#echo sliceColumn(gtMtx,56)
#echo gtMtx.len

#echo findCellPairs(gtMtx = gtMtx)
#var nSnpsCells = newSeqWith(ncells,0)
#selectTemplateCell
echo selectTemplateCell(gtMtx = gtMtx, nPairs =3)