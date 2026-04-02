######################################
# wrapper for rule: qc_assembly_transrate
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: qc_assembly_transrate \n##\n")
f.close()

version = str(subprocess.Popen("conda list", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

trinity = snakemake.input.trinity if hasattr(snakemake.input,'trinity') and len(snakemake.input.trinity)>0 else []
spades  = snakemake.input.spades  if hasattr(snakemake.input,'spades') and len(snakemake.input.spades) >0 else []
megahit = snakemake.input.megahit if hasattr(snakemake.input,'megahit') and len(snakemake.input.megahit)>0 else []

command = "transrate"+ \
          " --assembly "+snakemake.input.merged+","+",".join(trinity+spades+megahit)+ \
          " --output "+snakemake.params.prefix+ \
          " --threads "+str(snakemake.threads)+ \
          " --left "+snakemake.input.r1+(" --right "+snakemake.input.r2 if hasattr(snakemake.input,'r2') else "")+ \
          (" --reference "+snakemake.params.ref if snakemake.params.ref != "" else "")+ \
          " >> "+snakemake.log.run+" 2>&1"
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)

command = "cat "+snakemake.output.summ+" | awk -F, '{{n=split($1,s,\"/\");$1=s[n];print $0}}' OFS='\\t' > temp.$$ && mv temp.$$ "+snakemake.output.summ+" >> "+snakemake.log.run+" 2>&1"
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)
