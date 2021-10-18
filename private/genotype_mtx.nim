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


proc getGtMtx(ibam:Bam, ivcf:VCF, barcodeTable:TableRef, outGtMtx:FileStream, outSnpAnnot:FileStream, 
              outTotalCountMtx:FileStream, outAltCountMtx:FileStream,maxTotalReads:int,
              minTotalReads:int,  mapq: int,minbsq:int, barcodeTag:string, minsnpdepth:int,minCellDp:int, maxCellDp:int,p0 = 0.3, p1 = 0.8):int = 
  var variantIndex = 1
  var writtenVarintIndex = 1
  var calt,cellDP,geno,cellIndex,nnsize = 0
  var delBarcodes:seq[string]
  var cellGtNodes: Table[string, GtNode]
  var cellGtNodeCopy: Table[string, GtNode]
  for rec in ivcf.query("*"): 
    var rec_alt = rec.ALT[0][0]
    var rec_ref = rec.REF[0]
    cellGtNodes = findGtNodes(rec=rec, variantIndex = variantIndex, ibam=ibam, maxTotalReads=maxTotalReads,
                                  minTotalReads=minTotalReads, mapq = mapq, barcodeTable = barcodeTable,
                                  minbsq=minbsq,barcodeTag=barcodeTag)
    cellGtNodeCopy = cellGtNodes
    for cellbarcode in cellGtNodeCopy.keys:
      cellDP = cellGtNodes[cellbarcode].alleles.len
      if cellDP < minCellDp:
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
  echo "wrote " & $(writtenVarintIndex) & " variants to matrices and snpAnnot file"
  for outFileStream in [outGtMtx,outTotalCountMtx,outAltCountMtx]:
    outFileStream.setPosition(49)
    outFileStream.write($(writtenVarintIndex) & " " & $barcodeTable.len & " " & $nnsize) 
  for outFileStream in[outGtMtx,outTotalCountMtx,outAltCountMtx, outSnpAnnot]:
    outFileStream.close()
  return 0

let vcf_file = "data/swapped/FVB_NJ.mgp.v5.snps.dbSNP142.homo.alt.modified.swapped.GT.chr1.vcf.gz"
let bam_file = "data/WC_CNV_42/WC_CNV_42.bam"
let barcode_file = "data/WC_CNV_42/WC_CNV_42_filteredBC.tsv"
var
  ibam:Bam
  ivcf:VCF
  outGtMtx, outSnpAnnot, outTotalCountMtx, outAltCountMtx:FileStream
if not open(ibam, bam_file, threads=1, index = true):
    quit "couldn't open input bam"
if not open(ivcf, vcf_file, threads=1):
    quit "couldn't open: vcf file"

var hf = hts.hts_open(cstring(barcode_file), "r")
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

let s_Chrs = @["chr1"]
var chrom = "chr1"
let out_prefix = "test_data/getMtx/"
echo "generating gtMtx"
try: 
  outGtMtx = openFileStream(out_prefix & chrom &  "_gtMtx.mtx", fmWrite)
  outSnpAnnot = openFileStream(out_prefix & chrom &  "_snpAnnot.txt", fmWrite)
  outTotalCountMtx = openFileStream(out_prefix & chrom &  "_totalMtx.mtx", fmWrite)
  outAltCountMtx = openFileStream(out_prefix & chrom &  "_AltMtx.mtx", fmWrite)
except:
  stderr.write getCurrentExceptionMsg()

## write headers to those mtx files and the first line place holder for total_row total_column total_entry
let sparseMatrixHeader = "%%MatrixMarket matrix coordinate integer general"

for outFileStream in [outGtMtx,outTotalCountMtx,outAltCountMtx]:
  outFileStream.writeLine(sparseMatrixHeader)
  outFileStream.writeLine(' '.repeat(50))

outSnpAnnot.writeLine(join(["POS","REF", "ALT"], sep = "\t"))

discard getGtMtx(ibam = ibam, ivcf = ivcf, barcodeTable = barcodeTable, outGtMtx = outGtMtx, outTotalCountMtx = outTotalCountMtx ,outAltCountMtx = outAltCountMtx,
                 outSnpAnnot = outSnpAnnot, maxTotalReads = 30,
                 minTotalReads = 5, mapq = 20, minbsq = 13, minCellDp = 2,maxCellDp =10,barcodeTag = "CB",minsnpdepth=1)


var outGtMtxFile, outTotalCountMtxFile,outAltCountMtxFile : string

## sort entries in mtx files 
for chrom in s_Chrs:
  outGtMtxFile = out_prefix & chrom & "_gtMtx.mtx"
  outTotalCountMtxFile = out_prefix & chrom &  "_totalMtx.mtx"
  outAltCountMtxFile = out_prefix & chrom &  "_AltMtx.mtx"
  for mtxFile in [outGtMtxFile, outTotalCountMtxFile,outAltCountMtxFile]:
    var imtx = readMtx(mtx_file = mtxFile)
    discard sortWriteMtx(imtx, mtx_file = mtxFile)

