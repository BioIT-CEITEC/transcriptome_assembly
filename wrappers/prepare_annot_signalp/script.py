###########################################
# wrapper for rule: prepare_annot_signalp
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_signalp \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = snakemake.params.tool+"-register "+snakemake.input.licence+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

version = str(subprocess.Popen(snakemake.params.tool+" -version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

# (time signalp -f short -T tmp -n signalp.gff -m signalp.fa -v -l signalp.log ../termites_merged.concat_assemblies.okay.concat_assemblies.tr.transdecoder_dir/longest_orfs.pep ) > signalp.short.out
if snakemake.params.tool == "signalp4":
    command = "$(which time) "+snakemake.params.tool+" -f short"+\
              " -T " + snakemake.params.tmp + \
              " -n " + snakemake.output.ref + \
              " -m " + snakemake.output.fa + \
              " -v -l " + snakemake.log.run + \
              " " + snakemake.input.pep + \
              " > " + snakemake.output.tab + \
              " 2>> " + snakemake.log.run
elif snakemake.params.tool == "signalp5":
    command = "$(which time) "+snakemake.params.tool+\
              " -format short"+\
              " -tmp " + snakemake.params.tmp + \
              " -gff3 " + \
              " -mature " + \
              " -plot png" + \
              " -org " + snakemake.params.organism + \
              " -prefix "+snakemake.params.prefix + \
              " -fasta " + snakemake.input.pep + \
              " -verbose" + \
              " -stdout" + \
              " > " + snakemake.output.tab + \
              " 2>> " + snakemake.log.run
else:
    raise Exception("Unsupported version of SignalP!")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.tool == "signalp5":
    command = "mv "+snakemake.params.prefix+"_mature.fasta "+snakemake.output.fa+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
