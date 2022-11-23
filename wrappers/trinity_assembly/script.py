######################################
# wrapper for rule: trinity_assembly
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

## NOTE: there was a necessary change inside script /mnt/ssd/ssd_1/snakemake/.snakemake/conda/798d548e/opt/trinity-2.8.5/util/misc/plot_strand_specificity_dist_by_quantile.Rscript
## line: c = cumsum(data$total_reads)
## into: c = cumsum(as.numeric(data$total_reads))

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: trinity_assembly \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# trinity_bin_dir = os.path.dirname(str(subprocess.Popen("which Trinity", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8'))

if snakemake.params.use_genome:
    input_files = " --genome_guided_bam "+snakemake.input.bam
    additional = " --genome_guided_max_intron "+str(snakemake.params.max_intron)+" --genome_guided_min_coverage "+str(snakemake.params.min_map_cov)
    if snakemake.params.jaccard_clip:
        additional += " --jaccard_clip"
else:
    input_files = " --seqType fq --left " + snakemake.input.r1 + " --right " + snakemake.input.r2
    additional = ""

strandness = "" if snakemake.params.stranded == "" else " --SS_lib_type "+snakemake.params.stranded

command = "$(which time) Trinity" + \
               input_files + \
               " --min_contig_length " + str(snakemake.params.min_contig_len) + \
               " --output " + os.path.dirname(snakemake.output.fasta) + \
               " --workdir " + os.path.dirname(snakemake.output.fasta) + \
               " --KMER_SIZE " + str(snakemake.wildcards.kmer) + \
               " --min_kmer_cov " + str(snakemake.params.min_kmer_cov) + \
               " --CPU " + str(snakemake.threads) + \
               " --max_memory " + str(snakemake.resources.mem) + "G" + \
               strandness + \
               additional + \
               " >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if not snakemake.params.use_genome:
    command = "mv " + os.path.dirname(snakemake.output.fasta) + "/Trinity.fasta " + snakemake.output.fasta + " >> " + snakemake.log.run + " 2>&1"
else:
    command = "mv " + os.path.dirname(snakemake.output.fasta) + "/Trinity-GG.fasta " + snakemake.output.fasta + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if not snakemake.params.use_genome:
    command = "mv " + os.path.dirname(snakemake.output.fasta) + "/Trinity.fasta.gene_trans_map " + snakemake.output.gtm + " >> " + snakemake.log.run + " 2>&1"
else:
    command = "mv " + os.path.dirname(snakemake.output.fasta) + "/Trinity-GG.fasta.gene_trans_map " + snakemake.output.gtm + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# cleaning of assembly data
if not snakemake.params.use_genome:
    command = "rm -rf " + os.path.dirname(snakemake.output.fasta) + "/{{read_partitions,chrysalis,scaffolding_entries.sam,jellyfish.kmers.fa,*.ok,.*.ok,*cmd*,*.bt2,both*,inchworm*,partitioned_reads*,"+snakemake.params.project_d+",insilico_read_normalization/tmp_normalized_reads}}"+" >> "+snakemake.log.run+" 2>&1"
else:
    command = "find " + os.path.dirname(snakemake.output.fasta)+"/*"+\
              " -not -name "+os.path.basename(snakemake.output.fasta)+\
              " -not -name "+os.path.basename(snakemake.output.gtm)+\
              " -delete >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
