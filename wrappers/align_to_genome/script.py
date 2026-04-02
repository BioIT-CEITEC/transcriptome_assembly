######################################
# wrapper for rule: align_to_genome
######################################
import os
import subprocess
import glob
from snakemake.shell import shell

shell.executable("/bin/bash")

log_filename = str(snakemake.log.run)
f = open(log_filename, 'a+')
f.write("\n##\n## RULE: align_to_genome \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "mkdir -p "+snakemake.params.prefix+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

star_index_dir = snakemake.input.index[0].replace("/SAindex","")

f = open(log_filename, 'at')
strandness = False if snakemake.params.strandness == "unstr" else True
if strandness == True:
    extra_flags_star_motif = "" # For STAR outSAMstrandField
    # extra_flags_star_wig = " --outWigStrand Stranded" # For START bedGraph
    f.write("Running as Stranded experiment \n")
else:
    extra_flags_star_motif = " --outSAMstrandField intronMotif" # For STAR outSAMstrandField
    # extra_flags_star_wig = " --outWigStrand Unstranded" # For START bedGraph
    f.write("Running as Unstranded experiment \n")
f.close()


# if snakemake.params.paired == "SE":
#     STAR_parameters = " --chimSegmentMin 30"
# else:
#     STAR_parameters = " --peOverlapMMp 0.1 --chimOutJunctionFormat 1 --chimSegmentMin 12 --chimJunctionOverhangMin 12"

command = "$(which time) STAR --runMode alignReads --runThreadN " + str(snakemake.threads) + \
           " --genomeDir " + star_index_dir + \
           " --readFilesIn " + " ".join(snakemake.input.fastqs)  + \
           " --readFilesCommand zcat" + \
           " --sjdbOverhang " + str(snakemake.params.read_len) + \
           " --sjdbGTFfile " + str(snakemake.input.gtf) + \
           " --outFileNamePrefix " + snakemake.params.prefix + \
           " --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1" + \
           " --outFilterMismatchNmax " + str(snakemake.params.num_mismatch) + \
           " --outFilterMismatchNoverReadLmax 1.0" + \
           " --outFilterMismatchNoverLmax "+ str(snakemake.params.perc_mismatch) + \
           " --alignIntronMin 20 --alignIntronMax " + str(snakemake.params.max_intron) + \
           " --alignMatesGapMax " + str(snakemake.params.max_mate_dist) +\
           " --outFilterMatchNmin 0 --outFilterScoreMinOverLread " + str(snakemake.params.map_score)+\
           " --outFilterMatchNminOverLread " + str(snakemake.params.map_perc)+ \
           " --outSAMheaderHD @HD VN:1.4 SO:coordinate"+\
           " --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All" + \
           extra_flags_star_motif +" --sjdbScore 1 --twopassMode Basic " + \
           " --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate"+\
           " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)

command = "mv " + snakemake.params.prefix + "Aligned.sortedByCoord.out.bam " + snakemake.output.bam + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -rf " + snakemake.params.prefix + "*pass1" + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam + " " + snakemake.output.bai + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "(time STAR --runMode inputAlignmentsFromBAM" + \
#           " --inputBAMfile " + snakemake.output.bam + \
#           " --outWigType bedGraph" + extra_flags_star_wig + \
#           " --outFileNamePrefix " + snakemake.params.prefix+\
#           ") >> "+log_filename+" 2>&1 " # --outWigReferencesPrefix chr suitable for UCSC
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

