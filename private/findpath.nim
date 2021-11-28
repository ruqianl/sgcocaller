## implements the traceback function
import utils
import graph 
import streams
import tables
import math
import strutils


proc pathTrackBackIthSperm(currentSperm: var SpermViNodes,
                    ithSperm: int,
                    thetaRef: float,
                    thetaAlt: float,
                    cmPmb: float,
                    outFileVStateMtx: var FileStream,
                    viSegmentInfo: var FileStream,
                    posEnd: var int64,
                    inferProb: var float,
                    reverseProb: var float): int =
  var currentEm,prevEm:  array[stateRef..stateAlt, float]
  var posStart:int64
  var transitFlag = false
  var cSNP = 1
  var ithSNP: int
  var transProb: float
  var rightGap,leftGap: array[0..1, float]
  for i in 1..high(currentSperm.viNodeseq):
    var state =  currentSperm.viNodeseq[^i].state
    currentSperm.viNodeseq[^(i+1)].state = currentSperm.viNodeseq[^i].pathState[state]
    prevEm = getEmission(thetaRef=thetaRef,thetaAlt=thetaAlt,
                          cRef=currentSperm.viNodeseq[^(i+1)].cRef,
                          cAlt=currentSperm.viNodeseq[^(i+1)].cAlt)        
    ithSNP = currentSperm.spermSnpIndexLookUp[high(currentSperm.viNodeseq)-i+1]
    transProb = getTrans(currentSperm.viNodeseq[^(i+1)].pos,
                              currentSperm.viNodeseq[^(i)].pos,
                              cmPmb=cmPmb)
    if currentSperm.viNodeseq[^(i+1)].state == stateRef:
      outFileVStateMtx.writeLine($ithSNP & " " & $(ithSperm+1) & " 1")
      if state == stateRef:
        # not transitioning to a different state ie still in this segment of same state
        inferProb+=prevEm[stateRef]
        reverseProb+=prevEm[stateAlt]
        cSNP += 1
        transitFlag = false
        #posStart = currentSperm[^(i+1)].pos
      else: # there is transition to different state: ref(start) to alt(end) now output the segment info
        posStart = currentSperm.viNodeseq[^(i)].pos
        #leftGapSize = currentSperm[^(i)].pos - currentSperm[^(i+1)].pos
        leftGap=[math.ln(transProb),math.ln(1-transProb)]
        inferProb+=leftGap[0]
        reverseProb+=leftGap[1]
        viSegmentInfo.writeLine(join(["ithSperm" & $ithSperm, $posStart, $posEnd, $(inferProb-reverseProb) ,$cSNP, "2"],sep = " "))
        transitFlag = true
        cSNP = 1
        rightGap = leftGap
        inferProb = prevEm[stateRef]+rightGap[0]
        reverseProb = prevEm[stateAlt]+rightGap[1]
        posEnd = currentSperm.viNodeseq[^(i+1)].pos
    else:
      outFileVStateMtx.writeLine(join([$ithSNP, $(ithSperm+1), "2"],sep = " "))
      if state == stateAlt:
          # not transitioning to a different state ie still in this segment of same state
        inferProb += prevEm[stateAlt]
        reverseProb += prevEm[stateRef]
        cSNP += 1
        transitFlag = false
      else:
        ## state transit
        posStart = currentSperm.viNodeseq[^(i)].pos
        leftGap=[math.ln(transProb),math.ln(1-transProb)]
        inferProb += leftGap[0]
        reverseProb += leftGap[1]
        viSegmentInfo.writeLine(join(["ithSperm" & $ithSperm, $posStart, $posEnd, $(inferProb-reverseProb) ,$cSNP, "1"],sep = " "))
        transitFlag = true
        cSNP = 1
        rightGap = leftGap
        inferProb = prevEm[stateAlt]+rightGap[0]
        reverseProb = prevEm[stateRef]+rightGap[1]
        posEnd = currentSperm.viNodeseq[^(i+1)].pos
  ## traced to the start position of the chromosome for this cell 
  #leftGap = [0.0,0.0]
  if not transitFlag:
    ## the first SNP is included in the segment from traced from back
    posStart = currentSperm.viNodeseq[0].pos
    if currentSperm.viNodeseq[0].state == stateRef:
      viSegmentInfo.writeLine(join(["ithSperm" & $ithSperm, $posStart, $posEnd, $(inferProb-reverseProb) ,$cSNP, "1"],sep = " "))
    else:
      viSegmentInfo.writeLine(join(["ithSperm" & $ithSperm, $posStart, $posEnd, $(inferProb-reverseProb) ,$cSNP, "2"],sep = " "))
  else:
    ## the first node has different state from the second, the first node has its own segment
    posStart = currentSperm.viNodeseq[0].pos
    posEnd = posStart
    cSNP = 1
    currentEm = getEmission(thetaRef=thetaRef,thetaAlt=thetaAlt,cRef=currentSperm.viNodeseq[0].cRef, cAlt=currentSperm.viNodeseq[0].cAlt) 
    transProb = getTrans(currentSperm.viNodeseq[0].pos, currentSperm.viNodeseq[1].pos, cmPmb=cmPmb)
    if currentSperm.viNodeseq[0].state == stateRef:
      inferProb = currentEm[stateRef]+math.ln(transProb)
      reverseProb = currentEm[stateAlt]+math.ln(1-transProb)
      viSegmentInfo.writeLine(join(["ithSperm" & $ithSperm, $posStart, $posEnd, $(inferProb-reverseProb) ,$cSNP, "1"],sep = " "))
    else:
      inferProb = currentEm[stateAlt]+math.ln(transProb)
      reverseProb = currentEm[stateRef]+math.ln(1-transProb)
      viSegmentInfo.writeLine(join(["ithSperm" & $ithSperm, $posStart, $posEnd, $(inferProb-reverseProb) ,$cSNP, "2"],sep = " ") )
  
  return 0
proc pathTrackBack*(scSpermSeq:  var SeqSpermViNodes,
                    thetaRef: float,
                    thetaAlt: float,
                    cmPmb: float,
                    outFileVStateMtx: var FileStream,
                    viSegmentInfo: var FileStream): int = 

  var posEnd: int64
  var inferProb,reverseProb = 0.0
  var currentEm: array[stateRef..stateAlt, float]
  var lastNode: ViNode
  var spermVNseq: SpermViNodes
  var ithSNP: int

  for ithSperm in 0..(scSpermSeq.len-1):
    ## rightGap,leftGap = [0.0,0.0]
    spermVNseq = scSpermSeq[ithSperm]
    if spermVNseq.viNodeseq.len==0: continue
    lastNode = spermVNseq.viNodeseq[high(spermVNseq.viNodeseq)]
    currentEm = getEmission(thetaRef=thetaRef,thetaAlt=thetaAlt,cRef=lastNode.cRef, cAlt=lastNode.cAlt)
    if lastNode.pathScore[stateRef] > lastNode.pathScore[stateAlt]:
      scSpermSeq[ithSperm].viNodeseq[high(scSpermSeq[ithSperm].viNodeseq)].state = stateRef
      ithSNP = scSpermSeq[ithSperm].spermSnpIndexLookUp[high(scSpermSeq[ithSperm].viNodeseq)+1]
      outFileVStateMtx.writeLine($ithSNP & " " & $(ithSperm+1) & " 1")
      inferProb  = currentEm[stateRef]
      reverseProb = currentEm[stateAlt]
    else:
      scSpermSeq[ithSperm].viNodeseq[high(scSpermSeq[ithSperm].viNodeseq)].state = stateAlt
      ithSNP = scSpermSeq[ithSperm].spermSnpIndexLookUp[high(scSpermSeq[ithSperm].viNodeseq)+1]
      outFileVStateMtx.writeLine($ithSNP & " " & $(ithSperm+1) & " 2")
      inferProb  = currentEm[stateAlt]
      reverseProb = currentEm[stateRef]
    posEnd = lastNode.pos
    ## call pathTrackBack will write the most probably hidden state seq and viterbi segment info to the  relevant File Streams.
    discard pathTrackBackIthSperm(currentSperm = scSpermSeq[ithSperm], ithSperm = ithSperm,  thetaRef = thetaRef, 
                          thetaAlt = thetaAlt, cmPmb = cmPmb,outFileVStateMtx = outFileVStateMtx,
                          viSegmentInfo = viSegmentInfo, posEnd = posEnd,  inferProb = inferProb,reverseProb = reverseProb)                    

