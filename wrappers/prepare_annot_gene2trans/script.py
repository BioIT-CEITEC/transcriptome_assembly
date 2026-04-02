###########################################
# wrapper for rule: prepare_annot_gene2trans
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_gene2trans \n##\n")
f.close()


command = "(time grep '^>' " + snakemake.input.okay + " | sed 's/>//g' | awk 'BEGIN {{ FS=\" \"; OFS=\"\\t\"}} {{print $1, $1}}' " + \
          ") > " + snakemake.output.gene2trans + " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
