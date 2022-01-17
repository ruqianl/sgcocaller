## write output VCF file
import hts
import streams
import strutils

let debug = true

proc readNextPhased(phasedAnnotFS:var FileStream): seq[string] = 
  var phaseRec = newSeq[string](4)
  var lineRec: string
  while not phasedAnnotFS.atEnd():
    lineRec = phasedAnnotFS.readLine()
    phaseRec = lineRec.splitWhitespace()
    if debug and phaseRec.len != 4:
      quit "len not 4 at " & lineRec & $phasedAnnotFS.atEnd()
    if phaseRec[3] == "-1": continue
    else: break
  return phaseRec

proc writePhaseToVCF*(ivcfFile:string, ovcfFile: string, phasedAnnotFile: string, threads:int, add_header_string = """##sgcocaller_v0.1=phase""", chrom = "nil"):int = 

  if debug: echo "ivcffile " & ivcfFile
  if debug: echo "ovcfFile " & ovcfFile
  if debug: echo "phasedAnnotFile " & phasedAnnotFile
  var ivcf: VCF
  var ovcf: VCF
  var v_off = 0
  var phasedAnnotFS: FileStream
  var phaseRec = newSeq[string](4)
  var phasePos: int
  var phasePhase: int
  if not open(ivcf, ivcfFile, threads=threads):
      quit "couldn't open input vcf file"
  if not open(ovcf, ovcfFile, threads=threads, mode ="w"):
    quit "couldn't open output vcf file"
  try: 
    phasedAnnotFS = openFileStream(phasedAnnotFile, fmRead)
  except:
    stderr.write getCurrentExceptionMsg()
  ovcf.copy_header(ivcf.header)
  discard ovcf.header.add_string(add_header_string)
  discard ovcf.write_header()
  discard phasedAnnotFS.readLine()
  var gt_string:seq[int32]
  phaseRec = readNextPhased(phasedAnnotFS)
  if debug: echo phasedAnnotFS.atEnd()
  phasePos = parseInt(phaseRec[0])
  phasePhase = parseInt(phaseRec[3])
  var queryChrom = chrom
  if queryChrom == "nil":
    queryChrom = "*"
  if debug: echo "wrote headers, and now starting write phased VCF"
  for v in ivcf.query(queryChrom):
    if v.POS.cint != phasePos: continue
    else:
      var f = v.format()
#      if debug: echo "wrote v in ivcf and set GT " & $v.POS.cint
      if phasePhase == 0:
        ## 0|1
        #if debug: echo "block hap0 and this h0 is 0 0|1"
        gt_string = @[int32(2),int32(5)]
      elif phasePhase == 1:
        ## 1|0
        #if debug: echo "block hap0 and this h0 is 1 1|0"
        gt_string = @[int32(4),int32(3)]
      else:
        if phasedAnnotFS.atEnd(): break 
      if f.set("GT",gt_string) != Status.OK:
        quit "set GT failed"
      if not ovcf.write_variant(v) :
        quit "write vcf failed for " & $voff
      if phasedAnnotFS.atEnd(): break 
      phaseRec = readNextPhased(phasedAnnotFS)
      # if debug:
      #   if phaseRec.len != 4: 
      #     echo " phaseRec.len " &  $phaseRec.len
      #     echo " phaseRec" & $phaseRec
      phasePos = parseInt(phaseRec[0])
      phasePhase = parseInt(phaseRec[3])
  phasedAnnotFS.close()
  ovcf.close()
  ivcf.close()
  
  return 0

