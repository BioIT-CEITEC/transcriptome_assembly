######################################
# wrapper for rule: merge_assemblies
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: merge_assemblies \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

compact_fasta = os.path.abspath(snakemake.params.base + ".compact.fa")

command = "tar xf "+snakemake.input.evigene+" -C "+snakemake.params.tmpd+" >> "+ snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# cat $SCRATCH/$INPUT_FILE* > $SCRATCH/all_assemblies/All_assemblies.fasta
#
# #Rename fasta headers
# cd $SCRATCH/all_assemblies/
# perl -ane 'if(/\>/){$a++;print ">Locus_$a\n"}else{print;}' All_assemblies.fasta > All_renamed.fasta
#
# #Reformat fasta headers, then run the tr2aacds pipeline
# trformat.pl -output All_assemblies.tr -input All_renamed.fasta
# rm All_assemblies.fasta All_renamed.fasta
# tr2aacds.pl -mrnaseq All_assemblies.tr -NCPU=$THREADS 1>tr2aacds.log 2>tr2aacds.err

command = "rm -f " + compact_fasta + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

for inp in snakemake.input.assemblies:
  with open(snakemake.log.run, 'at') as f: f.write("# working on "+inp+"\n")
  if "megahit/intermediate_contigs" in inp:
    short_name = 'megahit_'+os.path.basename(inp).split('.')[0].replace('k','')
  else:
    short_name = os.path.basename(inp).split('.')[1]
  with open(snakemake.log.run, 'at') as f: f.write("# short name: "+short_name+"\n")
  # print("# short name: "+short_name)
  command = "cat " + inp + "|sed 's/^>/>"+short_name+"./'" \
            " >> " + compact_fasta + " 2>> " + snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

# command = "$(which time) --verbose cat " + " ".join(snakemake.input.trinity) + " " + snakemake.input.spades + \
#           " > " + snakemake.params.base + ".concat_assemblies.fa 2>> " + snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# command = "$(which time) --verbose "+snakemake.params.non_conda_tools_dir + "/evigene/scripts/rnaseq/trformat.pl" + \
#             " -output " + snakemake.params.base + ".compact.fa" + \
#             " -input " + snakemake.params.base + ".concat_assemblies.fa" + \
#             " >> " + snakemake.log.run + " 2>&1 "
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

command = "cd " + snakemake.params.base_dir + ";" + \
            "$(which time) " + os.path.abspath(snakemake.params.tmpd+"/evigene/scripts/prot/tr2aacds.pl") + \
            " -mrnaseq " + compact_fasta + \
            " -NCPU=" + str(snakemake.threads) + \
            " -MAXMEM="+str(snakemake.resources.mem)+"000" + \
            " -logfile -tidyup" + \
            " >> " + os.path.abspath(snakemake.log.run) + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.trclass+" "+snakemake.output.trclass+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -rf " + compact_fasta + " "+snakemake.params.base + ".compact{{,nrcd1_blsplit,_split,nrcd1_db.perf}} "+snakemake.params.base_dir+"/tmpfiles >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
