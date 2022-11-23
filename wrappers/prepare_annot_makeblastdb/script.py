###########################################
# wrapper for rule: prepare_annot_makeblastdb
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_makeblastdb \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

db_name = snakemake.input.fa[:-5]+"diamond_db"

# # Unpack, and build the downloaded Uniprot/Swissprot database
# command = "$(which time) --verbose unpigz -p " + str(snakemake.threads) + " -c " + snakemake.input.uniprot_db + \
#           " > " + snakemake.params.db_input + " 2>> " + snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

command = "$(which time) --verbose diamond makedb --in " + snakemake.input.fa+\
          " --db " + db_name+\
          " -p " + str(snakemake.threads)+\
          " >> " + snakemake.log.run + " 2>&1 && touch "+snakemake.output.db_ok
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "zcat "+snakemake.input.fa+" | sed '/^>/ s/sp|\([^|]\+\)|[^ \t]\+/\\1/' > "+snakemake.params.simple_names+\
          " 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
