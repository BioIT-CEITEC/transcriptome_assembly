###########################################
# wrapper for rule: multiqc_on_eval_PE
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: multiqc_on_eval_PE \n##\n")
f.close()

conda = str(subprocess.Popen("conda list 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+conda+"\n")
f.close()

#wd = snakemake.params.working_dir
#strand = snakemake.params.stranded # [RF, FR] for paired-end; RF is "typical" Illumina second strand, FR is first strand; [F, R] for single-end
#extra_flags_expression="--SS_lib_type "+strand if strand == "RF" or strand == "FR" else ""
#gene2trans = snakemake.input.gene2trans
#salmon_genes_tpm = snakemake.output.gene_tpm
#salmon_isoforms_tpm = snakemake.output.iso_tpm

configs = ""
for conf in snakemake.params.multiqc_configs:
  if os.path.isfile(conf):
    configs+= " -c "+conf

command = "$(which time) --verbose multiqc -f -d"+\
          " "+configs+\
          " -n "+snakemake.output.html+\
          " -x 'tmp' -x 'rules' -x 'wrappers' -x '.snakemake' -x '.taskrunner'"+\
          " "+snakemake.params.path+\
          " >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
