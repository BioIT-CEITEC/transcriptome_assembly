######################################
# wrapper for rule: read_correction
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: read_correction \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()


command = "mkdir -p "+os.path.dirname(snakemake.output.reads[0])+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

input_reads = " -C <(zcat "+snakemake.input.reads[0]+")"
if snakemake.params.paired and len(snakemake.input.reads) == 2:
  input_reads+= " <(zcat "+snakemake.input.reads[1]+")"
  
# Jellyfish counts - the GitHub version of RCorrector doesn't need this step but Conda does...
command = "$(which time) jellyfish count"+\
               " -t " + str(snakemake.threads) + \
               " -m " + str(snakemake.params.kmer) + \
               " -s 100M" + \
               " -o " + str(snakemake.params.mer_counts) + \
               input_reads + \
               " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) jellyfish dump " + \
               " -o " + snakemake.params.mer_counts_dump + \
               " " + snakemake.params.mer_counts + \
               " >> " + snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) jellyfish histo -t " + str(snakemake.threads) + \
               " -o " + snakemake.params.mer_counts + ".histo" + \
               " " + snakemake.params.mer_counts + \
               " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

input_reads = " -p "+snakemake.input.reads[0]
if snakemake.params.paired and len(snakemake.input.reads) == 2:
  input_reads+= " "+snakemake.input.reads[1]
  
command = "$(which time) rcorrector -t " + str(snakemake.threads) + \
               " -k " + str(snakemake.params.kmer) + \
               " -wk " + str(snakemake.params.wk) + \
               " -maxcor " + str(snakemake.params.max_corr) + \
               " -maxcorK " + str(snakemake.params.max_corr_k) + \
               input_reads + \
               " -c " + snakemake.params.mer_counts_dump + \
               " -od " + os.path.dirname(snakemake.output.reads[0]) + \
               " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

inp1 = os.path.dirname(snakemake.output.reads[0])+"/"+os.path.basename(snakemake.input.reads[0].replace("fastq.gz","cor.fq.gz"))
out1 = snakemake.output.reads[0][:-3]

if snakemake.params.paired and len(snakemake.input.reads) == 2:
  inp2 = os.path.dirname(snakemake.output.reads[1])+"/"+os.path.basename(snakemake.input.reads[1].replace("fastq.gz","cor.fq.gz"))
  out2 = snakemake.output.reads[1][:-3]
  
  command = "(paste <(zcat "+inp1+ ") <(zcat "+inp2+") | paste - - - - |"+\
            " awk  '$1 !~ /unfixable_error/ && $2 !~ /unfixable_error/ {{print $1,$3,$5,$7 > file1;"+\
            " print $2,$4,$6,$8 > file2}}' FS='\\t' OFS='\\n' file1="+out1+" file2="+out2+") >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "$(which time) pigz -p "+str(snakemake.threads)+" -f "+out1+" "+out2+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "rm -rf "+snakemake.params.mer_counts_dump+" "+inp1+" "+inp2+" >> "+snakemake.log.run+" 2>&1" 
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
else:
  command = "(unpigz -p "+str(snakemake.threads)+" -c  "+inp1+ " | paste - - - - |"+\
            " awk  '$1 !~ /unfixable_error/ {{print $1,$2,$3,$4 > file1}}' FS='\\t' OFS='\\n' file1="+out1+\
            ") >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "$(which time) pigz -p "+str(snakemake.threads)+" -f "+out1+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "rm -rf "+snakemake.params.mer_counts_dump+" "+inp1+" >> "+snakemake.log.run+" 2>&1" 
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  

