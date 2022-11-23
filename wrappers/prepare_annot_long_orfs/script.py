###########################################
# wrapper for rule: prepare_annot_long_orfs
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_long_orfs \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

strand = snakemake.params.stranded # [RF, FR] for paired-end; RF is "typical" Illumina second strand, FR is first strand; [F, R] for single-end
extra_flags_transdecoder="-S" if strand == "RF" or strand == "FR" else ""
transdecoder_dir = os.path.dirname(snakemake.output.orfs_pep)

command = "mkdir -p " + snakemake.params.prefix + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# Transdecoder
command = "cd " + snakemake.params.prefix + ";" + \
          " $(which time) TransDecoder.LongOrfs " + extra_flags_transdecoder + \
          " -t " + os.path.abspath(snakemake.input.fa) + \
          " --gene_trans_map "+os.path.abspath(snakemake.input.gene2trans) + \
          " -O " + os.path.abspath(transdecoder_dir) + \
          " >> " + os.path.abspath(snakemake.log.run) + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -rf "+snakemake.params.prefix+"{{pipeliner*,transdecoder_dir.__checkpoints_longorfs}} >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

