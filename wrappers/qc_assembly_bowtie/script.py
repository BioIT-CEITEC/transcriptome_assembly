###########################################
# wrapper for rule: qc_assembly_bowtie
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: qc_assembly_bowtie \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

####################
### OLD STUFF 
# command = "$(which time) bowtie-build --threads "+str(snakemake.threads)+" " + snakemake.input.okay + " " + snakemake.params.prefix + ".bowtie_index >> " + snakemake.log.run + " 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
# 
# command = "($(which time) bowtie -q --phred33-quals -n 2 -e 99999999 -l 25 -I 1 -X 1000 -p " + str(snakemake.threads) + " " + \
#           " -a -m 200 --chunkmbs 128 -S " + snakemake.params.prefix + ".bowtie_index " + \
#           " -1 " + snakemake.input.r1 + " -2 " + snakemake.input.r2 + " | samtools view -@ " + str(snakemake.threads) + " " + \
#           " -b - -o " + snakemake.params.tmp_bam + \
#           ") >> " + snakemake.log.run + " 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
####################

trinity_bin_dir = os.path.dirname(str(subprocess.Popen("which Trinity", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8'))

command = "$(which time) " + trinity_bin_dir + "/../opt/trinity*/util/misc/run_bowtie2.pl"+\
          " --target " + snakemake.input.fa + \
          " --left " + snakemake.input.r1 + \
          " --right " + snakemake.input.r2 + \
          " --CPU " + str(snakemake.threads) + \
          " > "+snakemake.params.bow_bam_tmp+\
          " 2>> "+ snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools sort"+\
          " -@ " + str(snakemake.threads) + \
          " -m "+str(math.floor(100/snakemake.threads))+"G"+\
          " -o " + snakemake.params.bow_bam + \
          " "+snakemake.params.bow_bam_tmp + \
          " >> " + snakemake.log.run + " 2>&1 && rm "+snakemake.params.bow_bam_tmp+" >> "+snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools flagstat"+\
          " -@ " + str(snakemake.threads) + \
          " " + snakemake.params.bow_bam + \
          " > " + snakemake.output.flagstat + \
          " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools idxstats"+\
          " -@ " + str(snakemake.threads) + \
          " " + snakemake.params.bow_bam + \
          " > " + snakemake.output.idxstats + \
          " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools stats"+\
          " -@ " + str(snakemake.threads) + \
          " " + snakemake.params.bow_bam + \
          " > " + snakemake.output.stats + \
          " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) picard CollectInsertSizeMetrics"+\
          " I=" + snakemake.params.bow_bam + \
          " O=" + snakemake.output.picard_txt + \
          " H=" + snakemake.output.picard_pdf + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) " + trinity_bin_dir + "/../opt/trinity*/util/misc/examine_strand_specificity.pl"+\
          " " + snakemake.params.bow_bam + \
          " " + snakemake.output.ss_analysis.replace(".dat.vioplot.pdf","") + \
          " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -f " + snakemake.input.fa + ".{{[0-9],rev}}* >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
