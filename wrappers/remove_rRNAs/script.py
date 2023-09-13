######################################
# wrapper for rule: remove_rRNAs
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: remove_rRNAs \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if snakemake.params.stranded == "RF":
  strandness = "--nofw" 
elif snakemake.params.stranded == "FR":
  strandness = "--norc" 
else:
  strandness = ""
  
if snakemake.params.paired and len(snakemake.input.reads) == 2:
  input_reads = " -1 " + snakemake.input.reads[0] + " -2 " + snakemake.input.reads[1]
else:
  input_reads = " -U " + snakemake.input.reads[0]

# bowtie2 --quiet --very-sensitive-local --phred33  -x $1 -1 $2 -2 $3 --threads 12 --met-file ${4}_bowtie2_metrics.txt --al-conc-gz blacklist_paired_aligned_${4}.fq.gz --un-conc-gz blacklist_paired_unaligned_${4}.fq.gz  --al-gz blacklist_unpaired_aligned_${4}.fq.gz --un-gz blacklist_unpaired_unaligned_${4}.fq.gz
command = "$(which time) bowtie2 --very-sensitive-local " + snakemake.params.phred + \
               " " + strandness + \
               " --threads " + str(snakemake.threads) + \
               " --met-file " + snakemake.log.meta + \
               " --un-conc-gz " + snakemake.params.pu + \
               " --un-gz " + snakemake.params.uu + \
               " -x " + snakemake.params.bw2_idx + \
               input_reads + \
               " > /dev/null" + \
               " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "$(which time) --verbose samtools view -b " + snakemake.params.aln + " > " + snakemake.output.aln + " 2>> " + snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
# 
# command = "rm " + snakemake.params.aln + " 2>> " + snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

