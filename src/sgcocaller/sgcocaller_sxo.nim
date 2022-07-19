## sgcocaller sxo, using pre-generated mtx files for finding crossovers using a HMM model
from findPath import pathTrackBack
from graph  import addViNodeIthSperm,SpermViNodes
# SeqSpermViNodes
import tables
import hts
import utils
import sequtils
import streams
# bins of SNP indexes
import strutils
import os
import math

proc readCountMtxToSeq*(mtxFileStream:FileStream, countMtx:var seq[seq[int16]], by_cell = false):int = 
  var currentLineSeq:seq[int]
  var currentLine:seq[string]
  ## i j 1-based from file
  if by_cell:
    while not mtxFileStream.atEnd():
      currentLine = mtxFileStream.readLine().splitWhitespace()
      currentLineSeq = map(currentLine, proc(x: string): int = parseInt(x))
      countMtx[(currentLineSeq[0]-1)][(currentLineSeq[1]-1)] = int16(currentLineSeq[2])
  else:
    while not mtxFileStream.atEnd():
      currentLine = mtxFileStream.readLine().splitWhitespace()
      currentLineSeq = map(currentLine, proc(x: string): int = parseInt(x))
      countMtx[(currentLineSeq[1]-1)][(currentLineSeq[0]-1)] = int16(currentLineSeq[2])
  return 0
proc addViNodeSXO(barcodeTable: TableRef, alleleCountTable: Table[string,allele_expr], scSpermSeq: TableRef[int,SpermViNodes],
                  snpIndex: int,
                  thetaRef: float,
                  thetaAlt: float,
                  rec_pos: int,
                  nnsize : var int,
                  initProb: array[stateRef..stateAlt, float],
                  cmPmb: float,
                  outaltCountMtxFS:FileStream, outtotalCountMtxFS:FileStream): int = 
  for bc, ac in alleleCountTable.pairs:
    var ithSperm = barcodeTable[bc]
    nnsize.inc()
    outaltCountMtxFS.writeLine($(snpIndex+1) & " " & $(ithSperm+1) & " " & $ac.calt)
    outtotalCountMtxFS.writeLine($(snpIndex+1) & " " & $(ithSperm+1) & " " & $(ac.calt+ac.cref))
    var emissionArray = getEmission(thetaRef=thetaRef,thetaAlt=thetaAlt,cRef=ac.cref,cAlt=ac.cAlt)
    discard addViNodeIthSperm(scSpermSeq = scSpermSeq, cAlt = ac.calt, cRef = ac.cref, ithSperm = ithSperm, emissionArray = emissionArray, snpIndex = (snpIndex+1), initProb = initProb,rec_pos =rec_pos,cmPmb = cmPmb)
  return 0
## barcodeTable, cell barcode:cell index
proc sgcocallerSXO*(barcodeTable:TableRef, phase_dir:string, out_dir:string, thetaREF:float, thetaALT:float, cmPmb:float,
                    s_Chrs:seq[string], initProb: array[stateRef..stateAlt, float], phasedSnpAnnotFileName:string,
                    batchSize:int): int =
  var ncells = barcodeTable.len
  var nsnps, ithSperm:int
  var currentEntrySeq: seq[string]
  var totalCountMtxByCell: seq[seq[int16]]
  var altCountMtxByCell: seq[seq[int16]]
  var alleleCountTable: Table[string,allele_expr]
  ## iterate through each selected chromosome
  for chrom in s_Chrs:
    var nnsize = 0
    ## number of non zeros = number of elements in the mtx
    var phasedSnpannoFS,totalCountMtxFS,altCountMtxFS,outFileVStateMtx,outaltCountMtxFS,outSnpAnnotFS, outtotalCountMtxFS,viSegmentInfo:FileStream
    let sparseMatrixHeader = "%%MatrixMarket matrix coordinate integer general"
    phasedSnpannoFS =  openFileStream(phasedSnpAnnotFileName, fmRead)
    try:
      #_totalCount
      totalCountMtxFS = openFileStream(phase_dir & chrom & "_totalMtx.mtx", fmRead)
      altCountMtxFS = openFileStream(phase_dir & chrom & "_altMtx.mtx", fmRead)
      outaltCountMtxFS = openFileStream(out_dir & chrom & "_altCount.mtx", fmWrite)
      outtotalCountMtxFS = openFileStream(out_dir & chrom & "_totalCount.mtx", fmWrite)
      outSnpAnnotFS = openFileStream(out_dir & chrom & "_snpAnnot.txt", fmWrite)
      outFileVStateMtx = openFileStream(out_dir & chrom & "_vi.mtx", fmWrite)
      viSegmentInfo = openFileStream(out_dir & chrom & "_viSegInfo.txt", fmWrite)
    except:
      stderr.write getCurrentExceptionMsg()    
    ## write headers to those mtx files and the first line place holder for total_row total_column total_entry    
    for fs in [outFileVStateMtx,outaltCountMtxFS,outtotalCountMtxFS]:
      fs.writeLine(sparseMatrixHeader)
      fs.writeLine(' '.repeat(50))  
    outSnpAnnotFS.writeLine(join(["POS","REF", "ALT"], sep = "\t"))
    discard totalCountMtxFS.readLine() 
    discard altCountMtxFS.readLine()
    discard altCountMtxFS.readLine()
    #N, i,j
    currentEntrySeq = totalCountMtxFS.readLine().splitWhitespace()
    nsnps = parseInt(currentEntrySeq[0])
    ## gtMtx is cell by Snp format
    totalCountMtxByCell = newSeqWith(nsnps,newSeq[int16](ncells))
    altCountMtxByCell = newSeqWith(nsnps,newSeq[int16](ncells))
    discard readCountMtxToSeq(mtxFileStream = totalCountMtxFS, countMtx = totalCountMtxByCell, by_cell = true)
    discard readCountMtxToSeq(mtxFileStream = altCountMtxFS, countMtx = altCountMtxByCell, by_cell = true)
    var isnpIndex = -1 
#    var osnpIndex = 0
    ## do crossover calling batch by batch to avoid using too much RAM
    var all_sperm_index = (0..(barcodeTable.len-1)).toSeq()
    let num_sub_batches = int(floor(all_sperm_index.len/batchSize))
    var sub_batches = all_sperm_index.distribute(num_sub_batches,spread=false)
    var scSpermSeq = newTable[int,SpermViNodes]()
    for sperm_batch in sub_batches:
#      echo "sperm_batch.len : " & $sperm_batch
      isnpIndex = -1 
#      osnpIndex = 0
      phasedSnpannoFS.setPosition(0)
      discard phasedSnpannoFS.readLine()
      while not phasedSnpannoFS.atEnd():
        currentEntrySeq = phasedSnpannoFS.readLine().splitWhitespace()
        isnpIndex.inc()
        if currentEntrySeq[3] == "-1":        
          continue
        else:
          alleleCountTable = initTable[string,allele_expr]()
          for bc,ithSperm in barcodeTable:
            if not (ithSperm in sperm_batch): continue
            if totalCountMtxByCell[isnpIndex][ithSperm] == 0:
              ## not adding this vi node to the spermSeq
              continue
            if currentEntrySeq[3] == "1":
              alleleCountTable[bc] = allele_expr(cref:(altCountMtxByCell[isnpIndex][ithSperm]), calt: (totalCountMtxByCell[isnpIndex][ithSperm] - altCountMtxByCell[isnpIndex][ithSperm]))
            else:
              alleleCountTable[bc] = allele_expr(cref:(totalCountMtxByCell[isnpIndex][ithSperm] - altCountMtxByCell[isnpIndex][ithSperm]), calt: (altCountMtxByCell[isnpIndex][ithSperm]))     
          if alleleCountTable.len == 0:
            continue
          # if currentEntrySeq[3] == "1":
          #   outSnpAnnotFS.writeLine(join([currentEntrySeq[0], currentEntrySeq[2], currentEntrySeq[1]], sep="\t") )
          # else:
          #   outSnpAnnotFS.writeLine(join([currentEntrySeq[0], currentEntrySeq[1], currentEntrySeq[2]], sep="\t") )
          discard addViNodeSXO(barcodeTable = barcodeTable,  alleleCountTable = alleleCountTable,  scSpermSeq = scSpermSeq,
                            thetaRef = thetaRef, thetaAlt = thetaAlt, snpIndex = isnpIndex, rec_pos = parseInt(currentEntrySeq[0]), nnsize = nnsize,
                            initProb = initProb, cmPmb = cmPmb, outaltCountMtxFS = outaltCountMtxFS, outtotalCountMtxFS = outtotalCountMtxFS)   
      discard pathTrackBack(scSpermSeq = scSpermSeq, thetaRef = thetaRef, thetaAlt=thetaAlt, cmPmb = cmPmb, outFileVStateMtx = outFileVStateMtx,
                            viSegmentInfo = viSegmentInfo)
      scSpermSeq.clear()
    phasedSnpannoFS.setPosition(0)
    discard phasedSnpannoFS.readLine()
    while not phasedSnpannoFS.atEnd():
      currentEntrySeq = phasedSnpannoFS.readLine().splitWhitespace()
      if currentEntrySeq[3] == "1":
        outSnpAnnotFS.writeLine(join([currentEntrySeq[0], currentEntrySeq[2], currentEntrySeq[1]], sep="\t") )
      elif currentEntrySeq[3] == "0":
        outSnpAnnotFS.writeLine(join([currentEntrySeq[0], currentEntrySeq[1], currentEntrySeq[2]], sep="\t") )
      else:
        outSnpAnnotFS.writeLine(join([currentEntrySeq[0], "*", "*"], sep="\t") )
    for fs in [outFileVStateMtx,outaltCountMtxFS,outtotalCountMtxFS]:
      fs.setPosition(49)
      fs.write($(isnpIndex+1) & " " & $barcodeTable.len & " " & $nnsize) 
    for outFileStream in [phasedSnpannoFS,totalCountMtxFS,altCountMtxFS,outFileVStateMtx,outaltCountMtxFS,outSnpAnnotFS, outtotalCountMtxFS,viSegmentInfo]:
      outFileStream.close()
  return 0
