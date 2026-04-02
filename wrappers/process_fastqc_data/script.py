######################################
# wrapper for rule: process_fastqc_data
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: process_fastqc_data \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# Consider removing over-represented sequences from FastQC, they should be BLASTed against SILVA rRNA DB and removed in positive case.
# To extract them use command: unzip -p ../input_files/raw_fastq_qc/termites_merged_R1_fastqc.zip */fastqc_data.txt | sed -n '/>>Overrep/,/>>END/{//b;p}' > cleaned_fastq/termites_merged_R2.overrep.tsv
# and then to use python or grep to identify such reads and get rid of them with their mates

input_reads = " -C <(zcat "+snakemake.input.reads[0]+")"
if snakemake.params.paired and len(snakemake.input.reads) == 2:
  input_reads+= " <(zcat "+snakemake.input.reads[1]+")"

command = "$(which time) jellyfish count"+\
          " -t "+str(snakemake.threads)+\
          " -m "+str(snakemake.params.kmer)+\
          " -s 100M"+\
          " -o "+snakemake.params.t1+\
          " --text"+\
          " -L 100"+\
          input_reads+\
          " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

## There could be used a check for how big portion of all reads the over-represented sequences are and filtered if this portion is below 0.1 % or so
# with snakemake.params.t1.open() as f1:
#     t1 = pd.read_csv(f1, '\t')
#     thr = snakemake.params.othr
#     if len(t1.loc[t1.Percentage > thr]) > 0:

if os.path.isfile(snakemake.params.t1) and os.path.getsize(snakemake.params.t1) > 0:
    command = "unzip -p "+snakemake.input.arch+" */fastqc_data.txt 2>> "+snakemake.log.run+" | awk '/^Total Sequences/ {{print $3}}' 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## Getting number of reads/fragments using command: "+command+"\n")
    f.close()

    total_seq = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    f = open(snakemake.log.run, 'at')
    f.write("## Total reads/fragments: "+total_seq+"\n")
    f.close()

    command = "cat "+snakemake.params.t1+\
              "| tail -n +2 | sort -k2,2nr"+\
              "| awk '$2 >= "+str(int(total_seq)*float(snakemake.params.low_lim))+" {{print \">seq_\"NR\"\\t\"$2,$1}}' OFS='\\n'"+\
              " > "+snakemake.params.f1+" 2>> "+snakemake.log.run
    # command = "cat "+snakemake.params.t1+" | grep -v '^#' | awk '{{printf \">overrep_seq_%s\\n%s\\n\",NR,$1}}' FS='\\t' > "+snakemake.params.f1+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    if os.path.isfile(snakemake.params.f1) and os.path.getsize(snakemake.params.f1) > 0:
        command = "$(which time) blastn"+\
                  " -query "+snakemake.params.f1+\
                  " -task blastn"+\
                  " -db "+snakemake.params.db+\
                  " -out "+snakemake.params.b1+\
                  " -outfmt '6 qaccver qlen qstart qend saccver slen sstart send evalue pident qcovs length mismatch gapopen bitscore'"+\
                  " -num_threads "+str(snakemake.threads)+\
                  " -max_target_seqs 5 >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "$(which time) Rscript "+snakemake.params.rscript+\
                  " "+snakemake.params.b1+\
                  " "+snakemake.params.f1+\
                  " "+snakemake.output.tsv+\
                  " "+snakemake.params.qcov+\
                  " "+snakemake.params.mism+\
                  " >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    else:
        command = "touch "+snakemake.output.tsv+" >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
else:
    command = "touch "+snakemake.output.tsv+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

