###########################################
# wrapper for rule: prepare_annot_build
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_build \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "mkdir -p " + snakemake.params.prefix + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cd " + snakemake.params.prefix + ";" + \
          " $(which time) Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate >> " + os.path.abspath(snakemake.log.run) + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# Use an HMM search to identify PFAM domains in the longest ORFs
command = "$(which time) unpigz -p " + str(snakemake.threads) + " " + snakemake.params.pfam_db + \
          " >> " + snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) hmmpress -f " + snakemake.output.pfam_db + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

