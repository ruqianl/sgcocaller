# phase blocks from cell fmf genotypes

import blocks_utils
#import readfmf
import tables
import math
import sequtils
import hts

let n_anchor_snps = 9
let posterior_threshold = 0.99
let confidence_threshold = 0.95
let debug = false

proc match_phase(cell_geno:seq[int], block_h0:seq[int]): int =
  if cell_geno.len != block_h0.len:
    quit "the haplotypes lengths should be equal" & " " & $cell_geno.len & " " & $cell_geno.len
  var diff = map(toSeq(0..high(cell_geno)), proc (x:int):int = abs(cell_geno[x]-block_h0[x]))
  #echo "matching " & $cell_geno & " to " & $block_h0 
  var posterior_phase = cal_posterior(n_snps = cell_geno.len, n_mis = sum(diff))
  ## unphased 
  if cell_geno.len < n_anchor_snps:
    if debug: echo "returned -1 because too few SNPs in the block covered by this cell. " & $diff.len
    return -1
  if float(sum(diff)) > (n_anchor_snps+1)/2:
    if posterior_phase > posterior_threshold:
      if debug: echo "returned 1 with diffs to h0 : " & $sum(diff) & " of " & $diff.len & ": posterior probs = " & $posterior_phase
      return 1
    else:
      if debug: echo "returned -1 because posterior prob is not high enough with diffs to h0 : " & $sum(diff) & " of " & $diff.len & ": posterior probs = " & $posterior_phase
      return -1
  else:
    if posterior_phase > posterior_threshold:
      if debug: echo "returned 0 with diffs to h0 : " & $sum(diff) & " of " & $diff.len & ": posterior probs = " & $posterior_phase
      return 0
    else:
      if debug: echo "returned -1 because posterior prob is not high enough with diffs to h0 : " & $sum(diff) & " of " & $diff.len & ": posterior probs = " & $posterior_phase
      return -1 

proc block_phase_in_cells*(cellGenoTable:TableRef, block_pos: seq[int], block_h0: seq[int], left:bool, cellBlockPhaseTable: TableRef[string, seq[int]]): int =
  var n_covered_snps = 0
  var covered_snps_cell_geno: seq[int]
  var covered_snps_block_geno: seq[int]
  var block_phase = 0
  for barcode in cellGenoTable.keys():
    n_covered_snps = 0
    covered_snps_cell_geno = newSeq[int]()
    covered_snps_block_geno = newSeq[int]()
    if left:
      if debug: echo "Phasing left of the block for cell "  & $barcode  
      for pos_i in 0..high(block_pos):
        if cellGenoTable[barcode].hasKey($block_pos[pos_i]) and block_h0[pos_i] != -1:
          #echo "using " & $block_pos[pos_i] & "for cell " & barcode
          n_covered_snps.inc
          covered_snps_cell_geno.add(cellGenoTable[barcode][$block_pos[pos_i]])
          covered_snps_block_geno.add(block_h0[pos_i])
        if n_covered_snps > n_anchor_snps:
          break
    else: 
      if debug: echo "Phasing right of the block for cell "  & $barcode  
      #  high(block_pos)..0
      for pos_i in countdown((block_pos.len - 1),0):
        if cellGenoTable[barcode].hasKey($block_pos[pos_i]) and block_h0[pos_i] != -1:
          #echo "using " & $block_pos[pos_i] & "for cell " & barcode
          n_covered_snps.inc
          covered_snps_cell_geno.add(cellGenoTable[barcode][$block_pos[pos_i]])
          covered_snps_block_geno.add(block_h0[pos_i])
        if n_covered_snps > n_anchor_snps:
          break   
    block_phase = match_phase(covered_snps_cell_geno,covered_snps_block_geno)
    if not cellBlockPhaseTable.hasKey(barcode):
      cellBlockPhaseTable[barcode] = newSeq[int]()
    cellBlockPhaseTable[barcode].add(block_phase)     
proc find_contig_blocks*(blocks:seq[Block],cellBlockLeftPhaseTable: TableRef[string, seq[int]],cellBlockRightPhaseTable: TableRef[string, seq[int]] ):seq[int] = 
  var contig_blocks = newSeq[int]()
  var count_missing =  newSeqWith[blocks.len,0]
  var contig_block_allow_missing = int(ceil(float(cellBlockLeftPhaseTable.len) * 0.1))
  # for b in 0..high(blocks):
  var left_right_phase: seq[tuple[left:int,right:int]]
  for bc in cellBlockLeftPhaseTable.keys():
    left_right_phase = zip(cellBlockLeftPhaseTable[bc],cellBlockRightPhaseTable[bc])
    if debug: echo bc & "\t" & $left_right_phase
    for bl in 0..high(left_right_phase):
      if left_right_phase[bl][0] == -1 or left_right_phase[bl][1] == -1 :
        count_missing[bl].inc
  if debug: echo "block missing in X number of cells:" & $count_missing
  for i, c in count_missing:
    if c <= contig_block_allow_missing: contig_blocks.add(i)
  return contig_blocks

## find linkage of two blocks
## 
proc find_linkage*(block1: int, block2: int, cellBlockLeftPhaseTable: TableRef[string, seq[int]],cellBlockRightPhaseTable: TableRef[string, seq[int]]): Block_linkage= 
  if debug: echo "linking " & $block1 & "and" & $block2
  var left_right_phase: seq[tuple[left:int,right:int]]
  var linkString: string
  if not (block1 < block2):
    quit "block1 needs to be left of block2"
  var blockLinkage = Block_linkage(prevBlock :block1, nextBlock: block2, switch: -1)
  for bc in cellBlockLeftPhaseTable.keys():
    left_right_phase = zip(cellBlockLeftPhaseTable[bc],cellBlockRightPhaseTable[bc])
    if debug: echo $left_right_phase[block1].right & $left_right_phase[block2].left
    linkString = $left_right_phase[block1].right & $left_right_phase[block2].left
    if linkString in @["00","11"]:
      blockLinkage.linkageType_00.inc
    elif linkString in @["01","10"]:
      blockLinkage.linkageType_01.inc
  blockLinkage.switch_posterior_prob = cal_switch_posterior(error_rate = 0.1, (blockLinkage.linkageType_00 + blockLinkage.linkageType_01), blockLinkage.linkageType_01)
  blockLinkage.switch_confidence =  log10(blockLinkage.switch_posterior_prob) - log10(1-blockLinkage.switch_posterior_prob)
  return blockLinkage
  
# contig_blocks, the blocks to build contigs. Rest of the blocks are filled in if can
proc build_contig*(contig_blocks:seq[int], cellBlockLeftPhaseTable: TableRef[string, seq[int]],cellBlockRightPhaseTable: TableRef[string, seq[int]] ):seq[Block_linkage] = 
  var blockLinkageSeq: seq[Block_linkage]
  var blockLinkage:Block_linkage
  for b in 0..(contig_blocks.len-2):
    blockLinkage = find_linkage(block1=contig_blocks[b], block2=contig_blocks[b+1],cellBlockLeftPhaseTable,cellBlockRightPhaseTable)
    if blockLinkage.switch_confidence < 0.0:
      blockLinkage.switch = 0
    else: 
      blockLinkage.switch = 1
    blockLinkageSeq.add(blockLinkage)
  return blockLinkageSeq

proc infer_block_linkage_to_one_contig_block*(contig_block:int,to_infer_block:int,  
                                             cellBlockLeftPhaseTable: TableRef[string, seq[int]],
                                             cellBlockRightPhaseTable: TableRef[string, seq[int]],
                                             at_left: bool):int = 
  if at_left:
    if debug: echo "linking " & $to_infer_block & "to contig block " & $contig_block & " from left"
    var blockLinkage = Block_linkage(prevBlock : to_infer_block, nextBlock :contig_block, switch : -1) 
    blockLinkage = find_linkage(block1 = to_infer_block, block2 = contig_block, cellBlockLeftPhaseTable,cellBlockRightPhaseTable)
    if abs(blockLinkage.switch_confidence) < confidence_threshold:
      return -1
    elif blockLinkage.switch_confidence > 0:
      return 1
    else:
      return 0
  else:
    if debug: echo "linking " & $to_infer_block & "to contig block " & $contig_block & " from right"
    var blockLinkage = Block_linkage(prevBlock: contig_block, nextBlock : to_infer_block , switch: -1) 
    blockLinkage = find_linkage(block1 = contig_block, block2 = to_infer_block, cellBlockLeftPhaseTable, cellBlockRightPhaseTable)
    if abs(blockLinkage.switch_confidence) < confidence_threshold:
      return -1
    elif blockLinkage.switch_confidence > 0:
      return 1
    else:
      return 0

proc infer_non_contig_blocks*(contig_blocks:seq[int], linkedContigs:seq[Block_linkage], blocks:seq[Block],cellBlockLeftPhaseTable: TableRef[string, seq[int]],cellBlockRightPhaseTable: TableRef[string, seq[int]]):seq[int] = 
  ## time to gather the non_contig forming blocks and fill their phase:
  ## Anchor using the first block in the contig blocks
  var hap_sw = 0
  var hap_sw_left = 0
  var hap_sw_right = 0
  
  var blocks_phase = newSeqWith(blocks.len,-1)
  var left_contig_block:int
  var right_contig_block: int 
  blocks_phase[contig_blocks[0]] = 0
  for lc in linkedContigs:
    if lc.switch == 0:
      if debug: echo "blocks_phase[lc.nextBlock] " & $blocks_phase[lc.nextBlock]
      
      blocks_phase[lc.nextBlock] = blocks_phase[lc.prevBlock] 
    else:
      blocks_phase[lc.nextBlock] = (1 xor blocks_phase[lc.prevBlock])  
  if debug: echo "block contig haplotype : " & $blocks_phase
  for b_j in 0..(blocks.len-1):
    if not (b_j in contig_blocks):
      if b_j < contig_blocks[0]:
        hap_sw = infer_block_linkage_to_one_contig_block(contig_block = contig_blocks[0], to_infer_block = b_j, cellBlockLeftPhaseTable,cellBlockRightPhaseTable, at_left=true)
        if hap_sw == 1:
          blocks_phase[b_j] = (1 xor contig_blocks[0]) 
        elif hap_sw == 0:
          blocks_phase[b_j] = contig_blocks[0] 
      elif b_j > contig_blocks[contig_blocks.len-1]:
        hap_sw = infer_block_linkage_to_one_contig_block(contig_block = contig_blocks[contig_blocks.len-1], to_infer_block = b_j,cellBlockLeftPhaseTable,cellBlockRightPhaseTable,at_left=false)
        if hap_sw == 1:
          blocks_phase[b_j] = (1 xor contig_blocks[0]) 
        elif hap_sw == 0:
          blocks_phase[b_j] = contig_blocks[0] 
      else:
        ## in between two contig blocks
        left_contig_block = b_j
        while not (left_contig_block in contig_blocks):
          left_contig_block = left_contig_block - 1
        right_contig_block = b_j
        while not (right_contig_block in contig_blocks):
          right_contig_block = right_contig_block+1
        hap_sw_left = infer_block_linkage_to_one_contig_block(contig_block = left_contig_block, to_infer_block = b_j, cellBlockLeftPhaseTable,cellBlockRightPhaseTable, at_left=false)
        hap_sw_right = infer_block_linkage_to_one_contig_block(contig_block = right_contig_block, to_infer_block = b_j, cellBlockLeftPhaseTable,cellBlockRightPhaseTable, at_left=true)
        if blocks_phase[left_contig_block] != blocks_phase[right_contig_block]:
          if hap_sw_right == -1:
            if hap_sw_left == 1:
              blocks_phase[b_j] = blocks_phase[right_contig_block]
            elif hap_sw_left == 0:
              blocks_phase[b_j] = blocks_phase[left_contig_block]
          elif hap_sw_right == 0 and (hap_sw_left == 1 or hap_sw_left == -1):
              blocks_phase[b_j] = blocks_phase[right_contig_block]
          elif hap_sw_right == 1 and (hap_sw_left == 0 or hap_sw_left == -1):
              blocks_phase[b_j] = blocks_phase[left_contig_block]
    else: continue 
  return blocks_phase 

proc write_linked_blocks*(vcfFile:string, outFile: string, blocks:seq[Block], blocks_phase:seq[int], threads:int):int = 
  var ivcf: VCF
  var ovcf: VCF
  var v_off = 0
  var current_block = 0
  if not open(ivcf, vcfFile, threads=threads):
      quit "couldn't open input vcf file"
  if not open(ovcf, outFile, threads=threads, mode ="w"):
    quit "couldn't open output vcf file"
  ovcf.header = ivcf.header
  discard ovcf.header.add_string("""##sgcocaller_v0.1=phaseBlocks""")
  discard ovcf.write_header()
  var gt_string:seq[int32]
  var block_pos_i = 0
  for v in ivcf.query("*"):
    v_off.inc
    if blocks_phase[current_block] == -1:
      current_block.inc
      if current_block == blocks.len:
        break
      block_pos_i = 0
      continue
    if v_off != blocks[current_block].snpPos[block_pos_i]:
      continue
    else: ## equal 
      var f = v.format()
      if blocks_phase[current_block] == 0:
        # don't need to flip
        if blocks[current_block].h0[block_pos_i] == 0:
          ## 0|1
          if debug: echo "block hap0 and this h0 is 0 0|1"
          gt_string = @[int32(2),int32(5)]
        else:
          ## 1|0
          if debug: echo "block hap0 and this h0 is 1 1|0"
          gt_string = @[int32(4),int32(3)]
      elif blocks_phase[current_block] == 1:
        if blocks[current_block].h0[block_pos_i] == 0:
          ## 1|0
          if debug: echo "block hap1 and this h0 is 0 1|0"
          gt_string = @[int32(4),int32(3)]
        else:
          ## 0|1
          if debug: echo "block hap1 and this h0 is 1 0|1"
          gt_string = @[int32(2),int32(5)]
      if f.set("GT",gt_string) != Status.OK:
        quit "set GT failed"
      if not ovcf.write_variant(v) :
        quit "write vcf failed for " & $voff
      block_pos_i.inc
    if block_pos_i == blocks[current_block].snpPos.len:
      current_block.inc
      if current_block == blocks.len:
        break
      block_pos_i = 0
  ovcf.close()
  ivcf.close()

