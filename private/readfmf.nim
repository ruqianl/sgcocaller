## read fmf 
## barcode -> snpPos: geno
import tables
import streams
import strutils

proc readFmf*(frag_file:string,cellGenoTable:TableRef):int = 
  var fmf: FileStream
  var n_blocks,line_index,b_pos : int 
  var current_frag, barcode, b_geno: string
  var current_frag_seq: seq[string]
  try: 
      fmf = openFileStream(frag_file, fmRead)
  except:
      stderr.write getCurrentExceptionMsg()
  while not fmf.atEnd():
    current_frag = fmf.readLine()
    current_frag_seq = current_frag.splitWhitespace()
    n_blocks = parseInt(current_frag_seq[0])
    barcode = current_frag_seq[1].split('.')[0]
    for b in 1..n_blocks:
      b_pos = parseInt(current_frag_seq[(b)*2])
      b_geno = current_frag_seq[(b)*2+1]
      for i, g in b_geno:
        if cellGenoTable.hasKey(barcode):
          discard cellGenoTable[barcode].hasKeyOrPut($(b_pos+i), parseInt($g))
        else:
          var geno_table = newTable[string,int]()
          geno_table[$(b_pos+i)] = parseInt($g)
          cellGenoTable[barcode] = geno_table
    line_index.inc
  echo "total fmfs " & $line_index
  fmf.close()


# discard readFmf("data/WC_CNV_42/fmf_snpdepth1/WC_CNV_42.chr16.fmf",cellGenoTable)
# for k in cellGenoTable.keys():
#   echo k
#   break
# echo cellGenoTable.len

# for key,value in cellGenoTable["TGTATTCAGGACAGCT-1"].pairs():
#   echo $key & "\t" & $value
