######################################
# wrapper for rule: megahit_assembly
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: megahit_assembly \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "rm -rf " + snakemake.params.odir + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.paired:
    input_files = " -1 " + snakemake.input.r1[0] + " -2 " + snakemake.input.r1[1]
else:
    input_files = " -r " + snakemake.input.r1[0]

# $(which time) --verbose megahit -1 corrected_fastq/termites_merged_R1.correct.fastq.gz -2 corrected_fastq/termites_merged_R2.correct.fastq.gz  --k-list 21,29,39,59,79,99,119,141 -m 0.4 -t 20 -o assembly/termites_merged/megahit/megahit_def_k-list --out-prefix termites_merged --tmp-dir /mnt/ssd/ssd_1/tmp/ 2>&1 | tee assembly/termites_merged/megahit/megahit_def_k-list.log
command = "$(which time) megahit " + \
               input_files + \
               " --k-list "+snakemake.params.kmers+\
               " --tmp-dir " + snakemake.params.tmpd + \
               " -t " + str(snakemake.threads) + \
               " -m " + str(snakemake.params.mem) + \
               " -o " + snakemake.params.odir + \
               " --out-prefix "+snakemake.params.prefix+\
               " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + snakemake.params.fasta + " " + snakemake.output.fasta + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
