######################################
# wrapper for rule: mark_duplicates
######################################
import os
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

log_filename = str(snakemake.log.run)
f = open(log_filename, 'a+')
f.write("\n##\n## RULE: mark_duplicates \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.log.mtx)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


if snakemake.params.mark_duplicates == True:
    if snakemake.params.UMI == "no_umi" or snakemake.params.umi_usage == "no":
        # command = "export LD_BIND_NOW=1"
        # f = open(log_filename, 'at')
        # f.write("## COMMAND: "+command+"\n")
        # f.close()
        # shell(command)

        command = "$(which time) picard MarkDuplicates"+\
                  " INPUT="+snakemake.input.bam+\
                  " OUTPUT="+snakemake.output.bam+\
                  " METRICS_FILE="+snakemake.log.mtx+\
                  " REMOVE_DUPLICATES="+str(snakemake.params.rmDup)+\
                  " ASSUME_SORTED=true"+\
                  " PROGRAM_RECORD_ID=null"+\
                  " VALIDATION_STRINGENCY=LENIENT"+\
                  " -Xmx"+str(snakemake.resources.mem)+"g"+\
                  " -Djava.io.tmpdir="+snakemake.params.tmpd+\
                  " >> "+log_filename+" 2>&1"

    else:
        command = "$(which time) umi_tools dedup"+\
                  " -I " + snakemake.input.bam +\
                  " -S " + snakemake.output.bam +\
                  " --log " + snakemake.log.mtx +\
                  " --temp-dir "+snakemake.params.tmpd+\
                  " --output-stats "+snakemake.params.umi_stats_prefix+\
                  " --extract-umi-method=read_id"+\
                  " --umi-separator='_'"+\
                  " --method=directional"+\
                  " --edit-distance-threshold=0"+\
                  " --spliced-is-unique"+\
                  " --multimapping-detection-method=NH"+\
                  " >> "+log_filename+" 2>&1"
                  
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "$(which time) samtools index -@" + str(snakemake.threads) + " " + snakemake.output.bam + " >> " + log_filename + " 2>&1 "
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
        
    if snakemake.params.keep_not_markDups_bam == False:
        command = "rm -f " + snakemake.input.bam+" "+ snakemake.input.bam + ".bai >> "+log_filename+" 2>&1"
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)
else:

    command = "mv -T " + snakemake.input.bam + " " + snakemake.output.bam + " >> " + log_filename + " 2>&1 "
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    command = "mv -T " + snakemake.input.bai + " " + snakemake.output.bai + " >> " + log_filename + " 2>&1 "
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
