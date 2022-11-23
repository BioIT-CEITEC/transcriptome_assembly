###########################################
# wrapper for rule: prepare_annot_transd_pred
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_transd_pred \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

transdecoder_dir = snakemake.params.transdecoder_dir

# Use an HMM search to identify PFAM domains in uniprot_sprot blastp hits
command = "cd " + transdecoder_dir + "; " + \
          "$(which time) TransDecoder.Predict"+\
          " -t " + os.path.abspath(snakemake.input.fa) + \
          " --retain_pfam_hits " + os.path.abspath(snakemake.input.pfam_out) + \
          " --retain_blastp_hits " + os.path.abspath(snakemake.input.blastp_out) + \
          " --output_dir " + os.path.abspath(transdecoder_dir) + \
          " >> " + os.path.abspath(snakemake.log.run) + " 2>&1"
print("## "+command)
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

