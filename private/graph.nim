import tables
import utils
import streams
import hts
import math

type 
  # Nucleotide = enum
  #   A, C, T, G    
  SpermViNodes* = object
    viNodeseq*: seq[ViNode]
    ## look up from SNPIndex to current Sperm SNP index
    snpIndexLookUp*: Table[int,int]
    ## look up from current Sperm SNP index to SNPIndex  
    spermSnpIndexLookUp*: Table[int,int]
  SeqSpermViNodes* = seq[SpermViNodes]


proc addViNode*(barcodeTable: TableRef, 
               alleleCountTable: Table[string,allele_expr],
               scSpermSeq: var SeqSpermViNodes,
               outFileTotalCountMtx: var FileStream,
               outFileAltCountMtx: var FileStream,
               nnsize: var int,
               mindp: int,
               maxdp: int,
               snpIndex: int,
               thetaRef: float,
               thetaAlt: float,
               rec_pos: int,
               initProb: array,
               cmPmb: float): int = 
  for bc, ac in alleleCountTable.pairs:
    # ac.tostring(acs)
    ## mindp, maxdp, they are values per cell
    if (ac.cref+ac.calt) <= mindp or (ac.cref+ac.calt) >= maxdp: continue        
    var ithSperm = barcodeTable[bc]
    ## write to mtx Ref count
    outFileTotalCountMtx.writeLine($snpIndex & " " & $(ithSperm+1) & " " & $(ac.cRef+ac.cAlt))
    outFileAltCountMtx.writeLine($snpIndex & " " & $(ithSperm+1) & " " & $ac.cAlt)
    nnsize += 1
    var emissionArray = getEmission(thetaRef=thetaRef,thetaAlt=thetaAlt,cRef=ac.cRef,cAlt=ac.cAlt)
    if scSpermSeq[ithSperm].viNodeseq.len==0:
      var currentViNode = ViNode()
      currentViNode.pathScore[stateRef] = math.ln(initProb[stateRef])+emissionArray[stateRef]
      currentViNode.pathScore[stateAlt] = math.ln(initProb[stateAlt])+emissionArray[stateAlt]
      currentViNode.pathState[stateRef] = stateN
      currentViNode.pathState[stateAlt] = stateN
      currentViNode.state = stateN
      currentViNode.pos = int(rec_pos)
      currentViNode.cAlt = int(ac.cAlt)
      currentViNode.cRef = int(ac.cRef)
      
      scSpermSeq[ithSperm].viNodeseq.add(currentViNode)
      scSpermSeq[ithSperm].snpIndexLookUp[snpIndex] = scSpermSeq[ithSperm].viNodeseq.len
      scSpermSeq[ithSperm].spermSnpIndexLookUp[scSpermSeq[ithSperm].viNodeseq.len] = snpIndex
    else:
      let preVNode = scSpermSeq[ithSperm].viNodeseq[high(scSpermSeq[ithSperm].viNodeseq)]
      var ltransProb = math.ln(getTrans(preVNode.pos,int(rec_pos),cmPmb=cmPmb))
      var lnoTransProb = math.ln(1-getTrans(preVNode.pos,int(rec_pos),cmPmb=cmPmb))
      var currentViNode = ViNode()
        # ref/alt -> ref
      var refTref = preVNode.pathScore[stateRef] + lnoTransProb
      var altTref = preVNode.pathScore[stateAlt] + ltransProb
      if refTref > altTref:
        currentViNode.pathScore[stateRef] = refTref
        currentViNode.pathState[stateRef] = stateRef
      else:
        currentViNode.pathScore[stateRef] = altTref
        currentViNode.pathState[stateRef] = stateAlt
        # ref/alt -> alt
      var refTalt = preVNode.pathScore[stateRef] + ltransProb
      var altTalt = preVNode.pathScore[stateAlt] + lnoTransProb
        
      if refTalt > altTalt:
        currentViNode.pathScore[stateAlt] = refTalt
        currentViNode.pathState[stateAlt] = stateRef
      else:
        currentViNode.pathScore[stateAlt] = altTalt
        currentViNode.pathState[stateAlt] = stateAlt
      currentViNode.pathScore[stateAlt] += emissionArray[stateAlt]
      currentViNode.pathScore[stateRef] += emissionArray[stateRef]
      currentViNode.cAlt = int(ac.cAlt)
      currentViNode.cRef = int(ac.cRef)
      currentViNode.pos = int(rec_pos)
      scSpermSeq[ithSperm].viNodeseq.add(currentViNode)
      scSpermSeq[ithSperm].snpIndexLookUp[snpIndex] = scSpermSeq[ithSperm].viNodeseq.len
      scSpermSeq[ithSperm].spermSnpIndexLookUp[scSpermSeq[ithSperm].viNodeseq.len] = snpIndex
  return 0
# The size line  
