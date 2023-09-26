######################################
# wrapper for rule: qc_assembly_rnaquast
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: qc_assembly_rnaquast \n##\n")
f.close()

version = str(subprocess.Popen("conda list", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# This fix needs to be here until the release of rnaQUAST v2.2.4 or later will be available using conda
command = "export FILE=$(find $CONDA_PREFIX -type f -name UtilsPictures.py 2>> "+snakemake.log.run+") && echo \"Checking file $FILE if missing 'import math'\" >> "+snakemake.log.run+" && [[ ! $(grep -F \"import math\" $FILE) ]] && sed -i '2 a import math' $FILE >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

genemark_dir = os.path.abspath(os.path.join(snakemake.params.tmpd,"GeneMarkS-T"))

command = "tar xf "+snakemake.input.genemark+" -C "+snakemake.params.tmpd+" >> "+ snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

strandness = " -ss" if snakemake.params.stranded != "" else ""
command = "export PATH="+genemark_dir+":$PATH && $(which time) --verbose rnaQUAST.py "+ \
          " -c "+" ".join(snakemake.input.assemblies)+ \
          " -o "+snakemake.params.odir+ \
          " -t "+str(snakemake.threads)+ \
          strandness+ \
          " --gene_mark"+ \
          " >> " + snakemake.log.run + " 2>&1"
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)

# with open(snakemake.log.run, 'at')
#   f.write("## INFO: writing header for MultiQC into "+snakemake.params.full_rep+"\n")
# with open(snakemake.params.full_rep, 'wt') as r:
#   r.write("""# id: rnaquast_table
# # section_name: 'rnaQUAST summary table'
# # description: 'This table contains merged data from short_report.txt and comparison_output/basic_metrics.txt of rnaQUAST tool.'
# # plot_type: table
# # format: tsv
# # pconfig:
# #  title: 'rnaQUAST statistics'""")

# cat short_report.tsv <(less comparison_output/basic_metrics.txt | grep -vP '^( |$)' | sed -e 's/ [ ]\+/\t/g' -e 's/\t$//'|tail -4) > full_report.tsv
command = "cat "+snakemake.params.odir+"short_report.tsv <(less "+snakemake.params.odir+"comparison_output/basic_metrics.txt"+\
          " | grep -vP '^( |$)' | sed -e 's/ [ ]\+/\\t/g' -e 's/\\t$//'|tail -4)"+\
          " > "+snakemake.params.odir+"full_report.tsv 2>> "+snakemake.log.run
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)

# Rscript -e  'data.table::fwrite(t(data.table::fread("full_report.tsv", sep="\\t", header=F)),"",sep="\\t",quote=F,col.names=F,row.names=F)' > full_report_t.tsv
command = 'Rscript -e  \'data.table::fwrite(t(data.table::fread("'+snakemake.params.odir+'full_report.tsv", sep="\\\\t", header=F))'+\
          ',"'+snakemake.output.full_rep+'",sep="\\\\t",quote=F,col.names=F,row.names=F)\' >> '+snakemake.log.run+" 2>&1"
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)

command = "cp "+snakemake.params.odir+"comparison_output/Nx.png "+snakemake.params.odir+"comparison_output/Nx_all.png >> "+snakemake.log.run+" 2>&1"
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)

command = "cp "+snakemake.params.odir+"comparison_output/transcript_length.png "+snakemake.params.odir+"comparison_output/transcript_length_all.png >> "+snakemake.log.run+" 2>&1"
with open(snakemake.log.run, 'at') as f:
  f.write("## COMMAND: "+command+"\n")
shell(command)
