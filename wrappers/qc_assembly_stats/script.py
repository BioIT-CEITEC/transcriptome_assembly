###########################################
# wrapper for rule: qc_assembly_stats
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: qc_assembly_stats \n##\n")
f.close()

version = str(subprocess.Popen("conda list", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

### General Trinity stats
command = "$(which time) TrinityStats.pl"+\
          " " + snakemake.input.fa + \
          " > " + snakemake.output.trin_stats + \
          " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

