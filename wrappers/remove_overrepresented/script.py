######################################
# wrapper for rule: remove_overrepresented
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: remove_overrepresented \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

out_R1 = snakemake.output.reads[0][:-3]
if snakemake.params.paired and len(snakemake.input.reads) == 2:
  out_R2 = snakemake.output.reads[1][:-3]

  command = "(paste -d '\\t' <(zcat "+snakemake.input.reads[0]+"|paste - - - - ) <(zcat "+snakemake.input.reads[1]+"|paste - - - - )"+\
            "|$(which time) grep -vFf "+snakemake.input.tsv+\
            "|awk '{{print $1,$2,$3,$4 > out1; print $5,$6,$7,$8 > out2 }}' FS='\\t' OFS='\\n'"+\
            " out1="+out_R1+" out2="+out_R2+") >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "$(which time) pigz -p "+str(snakemake.threads)+" -f "+out_R1+" "+out_R2+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
else:
  command = "(zcat "+snakemake.input.reads[0]+"|paste - - - - "+\
            "|$(which time) grep -vFf "+snakemake.input.tsv+\
            "|awk '{{print $1,$2,$3,$4 > out1}}' FS='\\t' OFS='\\n'"+\
            " out1="+out_R1+") >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "$(which time) pigz -p "+str(snakemake.threads)+" -f "+out_R1+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

