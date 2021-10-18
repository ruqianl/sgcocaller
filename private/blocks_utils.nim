import streams
import strutils
import math

let block_sep = "********"
let split_threshold = 30.0
let mis_threshold = 30.0
#let version = "v0.3.3"

type
  Block* = ref object
    # the positions of SNP (index position from provided VCF)
    snpPos*: seq[int]
    # number of SNPs it spans
    length*: int
    # number of SNPs phased
    phased*: int 
    # number of base pairs it spans
    # span 
    span*: int
    # number of fragments
    fragments*: int
    # haplotype type 0, its bit-wise complementary is h1
    h0*: seq[int] 
    switch_error*: seq[float]
    mismatch_error*: seq[float]
    # 1 switch, 0 no switch comparing to the prev block, -1 means unlinked
    switch: int
    # int frags;
    # float SCORE, bestSCORE, lastSCORE;
    # int* slist; // ordered list of variants in this connected component
    # int lastvar; // index of the first and last variants in this connected component
    # int iters_since_improvement;


type Block_linkage* = ref object
  prevBlock*: int
  nextBlock*: int
  linkageType_00* :int
  linkageType_01*: int
  switch_posterior_prob*: float
  linkage_confidence*: float
  switch*: int


proc cal_posterior*(error_rate = 0.1, n_snps:int, n_mis:int): float = 
  var ll_h0 = (1-error_rate)^(n_snps-n_mis)*error_rate^n_mis
  var ll_h1 = (1-error_rate)^(n_mis)*error_rate^(n_snps-n_mis)
  var new_p = max(ll_h0/(ll_h0+ll_h1),ll_h1/(ll_h0+ll_h1))
  return new_p

proc cal_switch_posterior*(error_rate = 0.1, n_blocks:int, n_sw:int): float = 
  var ll_01 = error_rate^(n_blocks-n_sw)*(1-error_rate)^n_sw
  var ll_00 = error_rate^(n_sw)*(1-error_rate)^(n_blocks-n_sw)
  return ll_01/(ll_01+ll_00)

proc is_block_header(line_string:string): bool = 
  if line_string.splitWhitespace()[0] == "BLOCK:":
  # .len == 11:
    return true
  else:
    return false

## to honor the current format of HAPCUT2 the long blocks are not splitted
## We look at the 10th column: switch error score to split the blocks
proc count_blocks(block_file:string, by_header = false): int =
  var blocks = 1
  var bfs: FileStream
  var switch_error: float
  try: 
    bfs = openFileStream(block_file, fmRead)
  except:
    stderr.write getCurrentExceptionMsg()
  var current_line: string
  while not bfs.atEnd():
    current_line = bfs.readLine()
    if current_line.splitWhitespace()[0] == block_sep:
      continue
    if by_header:
      if current_line.is_block_header:
        blocks.inc
      continue
    else:
      if current_line.is_block_header:
        continue
      switch_error = parseFloat(current_line.splitWhitespace()[9])
      if switch_error < split_threshold :
        blocks.inc
  bfs.close()
  if by_header:
    return blocks - 1
  else:
    return blocks
  
proc add_record(line_record: string, current_block:Block):int = 
  var line_record_seq = line_record.splitWhitespace()
  let snpPos = parseInt(line_record_seq[0])
  let switch_error = parseFloat(line_record_seq[9])
  let mismatch_error = parseFloat(line_record_seq[10])
  var h0: int
  try:
    h0 =  parseInt(line_record_seq[1])
  except:
    h0 = -1
#  let discrete_prune = parseInt(line_record_seq[8])
  current_block.snpPos.add(snpPos)
  if mismatch_error < mis_threshold:
    h0 = -1
  current_block.h0.add(h0)
  current_block.switch_error.add(switch_error)
  current_block.mismatch_error.add(mismatch_error)

proc read_blocks*(block_file:string,  by_header = false):seq[Block] = 
  # start with block 0  
  var b_j = 0
  var bfs: FileStream
  var current_line: string
  var blockSeq = newSeq[Block]()
  var line_index = 0
  var switch_error: float
  try: 
    bfs = openFileStream(block_file, fmRead)
  except:
    stderr.write getCurrentExceptionMsg()
  ## First line in block file, not used, a block header line
  current_line = bfs.readLine()
  var currentBlock = Block()
  while not bfs.atEnd():
    current_line = bfs.readLine()
    if current_line.splitWhitespace()[0] == block_sep:
      continue
    if by_header:
      if current_line.is_block_header:
        blockSeq.add(currentBlock)
        currentBlock = Block()
        continue
      else:
        ## add current line of record to this block
        discard add_record(current_line,currentBlock)
    else:
      if current_line.is_block_header:
        continue
      else:
        switch_error = parseFloat(current_line.splitWhitespace()[9])
        if switch_error < split_threshold :
          blockSeq.add(currentBlock)
          currentBlock = Block()
        discard add_record(current_line,currentBlock)
  blockSeq.add(currentBlock)
  bfs.close()
  return blockSeq
