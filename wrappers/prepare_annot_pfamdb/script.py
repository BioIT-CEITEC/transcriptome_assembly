###########################################
# wrapper for rule: prepare_annot_pfamdb
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_pfamdb \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

trin_pfam_out = snakemake.output.pfam_out

command = "$(which time) hmmscan --cpu " + str(snakemake.threads) + \
          " --domtblout " + trin_pfam_out + \
          " -o " + trin_pfam_out + \
          " " + snakemake.input.pfam_db + \
          " " + snakemake.input.fa + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

