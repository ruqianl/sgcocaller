import tables
import strutils
import hts
import utils
import streams
# minSnpDepth, this SNP should at least be covered by two cells

proc findGtNodes*(rec:Variant, variantIndex: int, ibam:Bam, maxTotalReads:int, minTotalReads:int, mapq: int,
                  barcodeTable:OrderedTableRef,
                  minbsq:int,
                  barcodeTag:string, minCellDp:int, minSnpDepth:int ): Table[string,GtNode] =

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
    if aln.flag.unmapped or aln.mapping_quality.cint < mapq or aln.flag.dup: continue
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
     #   echo "adding ref geno to cell " & currentCB
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
  var calt = 0 
  var cellDP = 0
  var delBarcodes:seq[string]

  for cellbarcode in barcodedNodesTable.keys:
    cellDP = barcodedNodesTable[cellbarcode].alleles.len
    if cellDP < minCellDp:
      delBarcodes.add(cellbarcode)
      continue
    calt = barcodedNodesTable[cellbarcode].alleles.count("1")
    if calt/cellDP < 0.3:
      barcodedNodesTable[cellbarcode].genotype = 0
    elif calt/cellDP > 0.8:
      barcodedNodesTable[cellbarcode].genotype = 1
    else:
      delBarcodes.add(cellbarcode)
      continue
  for bc in delBarcodes:
    barcodedNodesTable.del(bc)
  if barcodedNodesTable.len < minSnpDepth:
    return initTable[string,GtNode]() 
  return barcodedNodesTable  

## write the .fmf output stream
proc writeToFMF*(twoNodesFrag: Fragment, fmf_out:FileStream): Fragment =
  if twoNodesFrag.node2.variantIndex - twoNodesFrag.node1.variantIndex == 1:
 #   echo "writting for cell " & twoNodesFrag.cellbarcode & " with 1 block"
    if $twoNodesFrag.node1.genotype == "" or $twoNodesFrag.node2.genotype == "":
      quit "Genotype unknown when writing out"
    fmf_out.writeLine(join(["1",twoNodesFrag.cellbarcode & "." & $twoNodesFrag.counter, $twoNodesFrag.node1.variantIndex, 
                            $twoNodesFrag.node1.genotype & $twoNodesFrag.node2.genotype,"~~"], sep = " "))
  else:
 #   echo "writting for cell " & twoNodesFrag.cellbarcode & " with 2 block"
    fmf_out.writeLine(join(["2",twoNodesFrag.cellbarcode & "." & $twoNodesFrag.counter, $twoNodesFrag.node1.variantIndex,
                           $twoNodesFrag.node1.genotype, $twoNodesFrag.node2.variantIndex,
                           $twoNodesFrag.node2.genotype,"~~"], sep = " "))
  twoNodesFrag.counter = twoNodesFrag.counter + 1 
  twoNodesFrag.node1 = twoNodesFrag.node2
  twoNodesFrag.node2 = nil
  return twoNodesFrag

proc updateCellFragment*(barcodedNodesTable:Table[string,GtNode],
                         cellFragmentsTable:TableRef[string,Fragment],
                         fmf_out:FileStream): TableRef[string,Fragment] = 

 # echo "updating for " & $barcodedNodesTable.len
  for cellbarcode in barcodedNodesTable.keys:
    # cellDP = barcodedNodesTable[cellbarcode].alleles.len
    # if cellDP < minCellDp: continue
    # calt = barcodedNodesTable[cellbarcode].alleles.count("1")
    # if calt/cellDP < 0.3:
    #   barcodedNodesTable[cellbarcode].genotype = 0
    # elif calt/cellDP > 0.8:
    #   barcodedNodesTable[cellbarcode].genotype = 1
    # else:
    #   continue
    if cellFragmentsTable.hasKey(cellbarcode):
      if (not cellFragmentsTable[cellbarcode].node2.isNil) or cellFragmentsTable[cellbarcode].node1.isNil :
        quit "Writing fmf went wrong, the second allele shoud be moved to first after getting one fragment"
      else:
        cellFragmentsTable[cellbarcode].node2 = barcodedNodesTable[cellbarcode]
        var newFrag = writeToFMF(cellFragmentsTable[cellbarcode],fmf_out)
     #   echo "writing out fmf for " & cellbarcode
        cellFragmentsTable[cellbarcode] = newFrag
    else:
      var newFragment = Fragment(node1: barcodedNodesTable[cellbarcode], cellbarcode: cellbarcode, counter:1)
      cellFragmentsTable[cellbarcode] = newFragment
  return cellFragmentsTable

