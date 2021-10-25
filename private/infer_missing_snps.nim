## infer missing snp in template cell by looking at the snp's link to other nearby snps in other cells

## author Ruqian Lyu
## rlyu@svi.edu.au
import sequtils
import math
import utils

let nAnchorSnps = 10

proc sliceColumn*(gtMtx:seq[seq[BinaryGeno]],colIndex:int): seq[BinaryGeno] = 
  map(gtMtx, proc(x:seq[BinaryGeno]):BinaryGeno = x[colIndex])

# cell_geno seq of length nAnchorSnps
# temp_geno seq of length nAnchorSnps
# return whether cell_geno matches with temp_geno in hap0 or hap1(bitwise complementary to temp_geno)
proc matchTemplateHap(cell_geno:seq[BinaryGeno], temp_geno:seq[BinaryGeno],error_rate = 0.1): seq[float] = 
  var nmatch = 0
  for i,g in cell_geno:
    if g == temp_geno[i]:
      nmatch += 1
  let mis = cell_geno.len - nmatch
  let ll_h0 = error_rate^mis*(1-error_rate)^(nmatch)
  let ll_h1 = error_rate^nmatch*(1-error_rate)^(mis)
  if ll_h0 > ll_h1:
    return @[0.0, ll_h0/(ll_h0 + ll_h1)]
  elif ll_h0 == ll_h1:
    return @[-1.0, 0.5]
  else:
    return @[1.0, ll_h1/(ll_h0 + ll_h1)]

proc inferGeno(type00:int, type10:int, error_rate = 0.1,posterior_thresh = 0.98): BinaryGeno = 
  let prob_0 = error_rate^type10*(1-error_rate)^(type00)
  let prob_1 = error_rate^type00*(1-error_rate)^(type10)
  let prob_0_p = prob_0/(prob_0 + prob_1)
  let prob_1_p = 1 - prob_0_p 
  if max(prob_0_p,prob_1_p) > posterior_thresh:
    if prob_0_p > prob_1_p:
      return gREF
    else:
      return gALT
  return gUnknown
# return full sequence of template geno by inferring missing SNPs's genotypes 
proc inferSnpGeno*(templateGeno: seq[BinaryGeno], gtMtx:seq[seq[BinaryGeno]], posterior_thresh = 0.99):seq[BinaryGeno] = 
  var fullGeno = templateGeno
  let nSnps = gtMtx[0].len
  var coexisPos:seq[int]
  var snpLD: seq[float]
  var type00, type10: int

  var snpiInCells = newSeq[BinaryGeno](gtMtx.len)
  var offset = 1
  for i,missingSnpi in templateGeno:
    if missingSnpi != gUnknown: continue
    snpiInCells = sliceColumn(gtMtx,i)
    type00 = 0
    type10 = 0
    for j,snpiGeno in snpiInCells:
      if snpiGeno == gUnknown: continue
      # find 10 coexisting positions for cell j with template geno seq
      offset = 1
      coexisPos = newSeq[int]()
      while(coexisPos.len < nAnchorSnps):
        if ((i+offset) < nSnps):
          if((gtMtx[j][i+offset] != gUnknown) and (fullGeno[i+offset] != gUnknown)):
            coexisPos.add(i+offset)
        if (i-offset) >= 0:
          if((gtMtx[j][i-offset] != gUnknown) and (fullGeno[i-offset] != gUnknown)):
            coexisPos.add(i-offset)
        offset += 1
      snpLD = matchTemplateHap(cell_geno = map(coexisPos,proc(x:int): BinaryGeno =  gtMtx[j][x]), 
                               temp_geno = map(coexisPos,proc(x:int): BinaryGeno =  fullGeno[x]) )
      if snpLD[1] > posterior_thresh:
        if snpLD[0] == 1.0:
          if snpiGeno == gREF: ## snpiGeno is 1 or 2 
            type10 += 1 
          else:
            type00 += 1
        else:
          if snpiGeno == gREF: 
            type00 += 1 
          else:
            type10 += 1
      else: continue
    # echo "type00 " & $type00
    # echo "type10 " & $type10
    fullGeno[i] = inferGeno(type00, type10)  
  return fullGeno

