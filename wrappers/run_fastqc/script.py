######################################
# wrapper for rule: run_fastqc
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: run_fastqc \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "$(which time) fastqc -o " + os.path.dirname(str(snakemake.input.fastq)) + \
               " " + snakemake.params.extra + \
               " -d " + snakemake.params.tmp + \
               " "+snakemake.input.fastq + \
               " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + snakemake.params.html + " " + snakemake.output.html + " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

