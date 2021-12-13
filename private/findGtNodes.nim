## Genrate GtNode per variant (SNP)

import tables
import hts
import utils

# minSnpDepth = 2, this SNP should at least be covered by two cells
proc findGtNodes*(rec:Variant, variantIndex: int, ibam:Bam, maxTotalReads:int, minTotalReads:int, mapq: int,
                  barcodeTable:TableRef,
                  minbsq:int,
                  barcodeTag:string ): Table[string,GtNode] =

  var barcodedNodesTable = initTable[string,GtNode]()
  var rec_alt:char
  var total_reads = 0
  rec_alt = rec.ALT[0][0]
  for aln in ibam.query(chrom = $rec.CHROM,start = rec.POS.cint-1, stop = rec.POS.cint):
    var cbt = tag[string](aln, barcodeTag)
    var currentCB: string
    var base_off: int
    var base: char
    if cbt.isNone:
      continue  
    else:
      currentCB = cbt.get      
    if aln.flag.unmapped or aln.mapping_quality.cint < mapq or aln.flag.dup:
#      echo "read skipped, bad mapping"
      continue
    ## skip unmapped, duplicated, mapq low reads or aln.flag.secondary or aln.flag.supplementary
    ## if not aln.flag.proper_pair: continue
    if not barcodeTable.hasKey(currentCB): 
      continue
    base_off = get_base_offset(position = rec.POS.cint, align = aln) 
    base = aln.base_at(base_off)
    if aln.base_quality_at(base_off).cint < minbsq: 
      continue
    total_reads+=1
    if barcodedNodesTable.hasKey(currentCB):
      if base == rec.REF[0]: 
    #     echo "adding ref geno to cell " & currentCB
        barcodedNodesTable[currentCB].alleles = add_allele(barcodedNodesTable[currentCB],alleleBinGeno = 0)
        continue
      if base == rec_alt: 
        barcodedNodesTable[currentCB].alleles = add_allele(barcodedNodesTable[currentCB],alleleBinGeno = 1)
        continue
    else:
    #  echo "creating new node for cell " & currentCB
      var newGtNode = GtNode(chrom: $rec.CHROM, variantIndex: variantIndex)
      if base == rec.REF[0]: 
        newGtNode.alleles = "0"
        barcodedNodesTable[currentCB] = newGtNode
        continue
      if base == rec_alt: 
        newGtNode.alleles = "1"
        barcodedNodesTable[currentCB] = newGtNode
        continue
  if total_reads > maxTotalReads or total_reads < minTotalReads:
    return initTable[string,GtNode]()
  return barcodedNodesTable  
