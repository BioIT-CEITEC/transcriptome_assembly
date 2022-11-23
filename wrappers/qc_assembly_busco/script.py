###########################################
# wrapper for rule: qc_assembly_busco
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: qc_assembly_busco \n##\n")
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

command = "$(which time) tar xzf " + snakemake.input.busco + " -C " + snakemake.params.prefix + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

lineage = str(snakemake.wildcards.lineage)

# BUSCO does not take output dir prefix with path
busco_dir = "busco_run"
command = "$(which time) busco" + \
          " -i " + snakemake.input.fa + " " + \
          " -o " + busco_dir + \
          " -l " + snakemake.params.prefix + lineage + \
          " -m tran" + \
          " -c " + str(snakemake.threads) + \
          " --out_path "+ snakemake.params.prefix + \
          " --download_path " + snakemake.params.tmpd + \
          " --force --tar --offline --quiet" + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

busco_path = snakemake.params.prefix + busco_dir + "/"
command = "rm -rf " + snakemake.params.prefix + lineage + \
          " "+busco_path + "short_summary*" + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

busco_path = busco_path + "run_" + lineage + "/"
command = "mv " + busco_path + "short_summary.txt" + \
          " "+snakemake.output.summ+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + busco_path + "full_table.tsv" + \
          " "+snakemake.output.ftab+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + busco_path + "missing_busco_list.tsv" + \
          " "+snakemake.output.miss+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
