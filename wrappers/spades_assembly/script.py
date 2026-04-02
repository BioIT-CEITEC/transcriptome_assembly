######################################
# wrapper for rule: spades_assembly
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: spades_assembly \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if snakemake.params.paired:
    input_files = " -1 " + snakemake.input.r1[0] + " -2 " + snakemake.input.r1[1]
    strandness  = " --ss rf " if snakemake.params.stranded.startswith("R") else " --ss fr "
else:
    input_files = " -s " + snakemake.input.r1[0]
    strandness  = ""

command = "export OMP_NUM_THREADS=" + str(snakemake.threads) + "; $(which time) spades.py " + \
               " --rna " + \
               input_files + \
               " --tmp-dir " + snakemake.params.tmpd + \
               " -t " + str(snakemake.threads) + \
               " -m " + str(snakemake.resources.mem) + \
               " -k " + str(snakemake.wildcards.kmer) + \
               strandness + \
               " -o " + snakemake.params.odir + \
               " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + snakemake.params.odir + "/transcripts.fasta " + snakemake.output.fasta + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

