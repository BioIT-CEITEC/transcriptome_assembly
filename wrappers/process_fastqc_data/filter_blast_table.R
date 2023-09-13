library(data.table)

args <- commandArgs(trailingOnly = T)

bt = args[1]
fa = args[2]
out = args[3]
qcov = as.numeric(args[4])
mism = as.numeric(args[5])

colnames = c("qaccver","qlen","qstart","qend","saccver","slen","sstart","send","evalue","pident","qcovs","length","mismatch","gapopen","bitscore")
bt = fread(bt, sep = "\t", col.names = colnames)

if(nrow(bt) == 0) {
  fwrite(data.table(), out, sep = '\t', row.names = F, col.names = F)

  } else {
  fa = fread(cmd=paste0("cat ",fa," | paste - - | sed 's/^>//' "), header = F, sep = '\t', col.names = c("qaccver","qcount","qseq"), key = "qaccver")
  
  # result of filtering input data from BLAST to avoid matches with low coverage or too many mismatches
  res = bt[length/qlen>=qcov & mismatch<=mism, .SD[order(evalue)][1], by=.(qaccver)]
  # write out over-represented sequences with strong evidence
  fwrite(fa[res, .(qseq)], out, sep = '\t', row.names = F, col.names = F)
}

