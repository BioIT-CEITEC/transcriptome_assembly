######################################
# wrapper for rule: merge_samples
######################################
import subprocess
import os
from snakemake.shell import shell

shell.executable("/bin/bash")

log_filename = str(snakemake.log.run)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_samples \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "($(which time) zcat "+" ".join(snakemake.input.reads)+" | pigz -c > "+snakemake.output.reads+" ) 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
