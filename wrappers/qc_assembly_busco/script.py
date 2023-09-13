###########################################
# wrapper for rule: qc_assembly_busco
###########################################
import os
import sys
import math
import subprocess
import uuid
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: qc_assembly_busco \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

workdir = os.path.join(snakemake.params.tmpd, "busco_"+uuid.uuid4().hex)
command = "mkdir -p " + snakemake.params.prefix + " " + workdir + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) tar xzf " + snakemake.input.busco + " -C " + workdir + " >> " + snakemake.log.run + " 2>&1"
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
          " -l " + os.path.join(workdir, lineage) + \
          " -m tran" + \
          " -c " + str(snakemake.threads) + \
          " --out_path "+ workdir + \
          " --download_path " + snakemake.params.tmpd + \
          " --force --tar --offline --quiet" + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

busco_path = os.path.join(workdir, busco_dir)
command = "rm -rf " + os.path.join(workdir, lineage) + \
          " "+busco_path + "/short_summary*" + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

busco_path = os.path.join(busco_path, "run_"+lineage)
command = "mv " + os.path.join(busco_path, "short_summary.txt") + \
          " "+snakemake.output.summ+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + os.path.join(busco_path, "short_summary.json") + \
          " "+snakemake.params.prefix+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + os.path.join(busco_path, "full_table.tsv") + \
          " "+snakemake.output.ftab+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + os.path.join(busco_path, "missing_busco_list.tsv") + \
          " "+snakemake.output.miss+" >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
