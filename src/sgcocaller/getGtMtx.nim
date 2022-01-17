## generate gtMtx, gtMtx_snpAnnot, totalCount.mtx altCount.mtx

## author Ruqian Lyu
## Date: 2021 Oct
## gtMtx in 1, or 2
import tables
import strutils 
import hts
import utils
import streams
import findGtNodes

## these outputs will be per chromosome

proc getGtMtx*(ibam:Bam, ivcf:VCF, barcodeTable:TableRef, outGtMtx:FileStream, outSnpAnnot:FileStream, 
              outTotalCountMtx:FileStream, outAltCountMtx:FileStream,maxTotalReads:int,
              minTotalReads:int, mapq: int,minbsq:int, barcodeTag:string, minsnpdepth:int,minCellDp:int, maxCellDp:int,p0 = 0.3, p1 = 0.7,chrom: string):int = 
  var variantIndex = 1
  var writtenVarintIndex = 1
  var calt,cellDP,geno,cellIndex,nnsize = 0
  var cellGtNodes,cellGtNodeCopy: Table[string, GtNode]
    
  ## write headers to those mtx files and the first line place holder for total_row total_column total_entry
  let sparseMatrixHeader = "%%MatrixMarket matrix coordinate integer general"

  for outFileStream in [outGtMtx,outTotalCountMtx,outAltCountMtx]:
    outFileStream.writeLine(sparseMatrixHeader)
    outFileStream.writeLine(' '.repeat(50))
  outSnpAnnot.writeLine(join(["POS","REF", "ALT"], sep = "\t"))
  for rec in ivcf.query(chrom): 
    var rec_alt = rec.ALT[0][0]
    var rec_ref = rec.REF[0]
    cellGtNodes = findGtNodes(rec=rec, variantIndex = variantIndex, ibam=ibam, maxTotalReads=maxTotalReads,
                                  minTotalReads=minTotalReads, mapq = mapq, barcodeTable = barcodeTable,
                                  minbsq=minbsq,barcodeTag=barcodeTag)
    cellGtNodeCopy = cellGtNodes
    for cellbarcode in cellGtNodeCopy.keys:
      cellDP = cellGtNodes[cellbarcode].alleles.len
      if cellDP < minCellDp or cellDP > maxCellDp:
        cellGtNodes.del(cellbarcode)
        continue
      calt = cellGtNodes[cellbarcode].alleles.count("1")
      if calt/cellDP < p0:
        cellGtNodes[cellbarcode].genotype = gREF
      elif calt/cellDP > p1:
        cellGtNodes[cellbarcode].genotype = gALT
      else:
        cellGtNodes.del(cellbarcode)
        continue
    if cellGtNodes.len < minSnpDepth:
      variantIndex.inc
      continue
    else:
      for cellbarcode in cellGtNodes.keys:
        cellDP = cellGtNodes[cellbarcode].alleles.len
        calt = cellGtNodes[cellbarcode].alleles.count("1")
        ## $cellGtNodes[cellbarcode].genotype BinaryGeno to "0", "1" or "-1" (not -1), write to file coding in 1 or 2
        geno = parseInt($cellGtNodes[cellbarcode].genotype) + 1
        # index from 1 for matrices
        cellIndex = barcodeTable[cellbarcode] + 1
        outGtMtx.writeLine(join([$writtenVarintIndex, $cellIndex, $geno],sep = " "))
        outTotalCountMtx.writeLine(join([$writtenVarintIndex, $cellIndex, $cellDP],sep = " "))
        outAltCountMtx.writeLine(join([$writtenVarintIndex, $cellIndex, $calt],sep = " "))
        nnsize.inc
      outSnpAnnot.writeLine(join([$rec.POS, $rec_ref, $rec_alt], sep="\t") )
      writtenVarintIndex.inc
    ## write to matrices 
          
  echo "processed " & $(variantIndex) & " variants"
  echo "wrote " & $(writtenVarintIndex-1) & " variants to matrices and snpAnnot file"
  for outFileStream in [outGtMtx,outTotalCountMtx,outAltCountMtx]:
    outFileStream.setPosition(49)
    outFileStream.write($(writtenVarintIndex-1) & " " & $barcodeTable.len & " " & $nnsize) 
  for outFileStream in[outGtMtx, outTotalCountMtx, outAltCountMtx, outSnpAnnot]:
    outFileStream.close()
  return 0
