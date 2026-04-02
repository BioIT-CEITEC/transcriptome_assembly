###########################################
# wrapper for rule: prepare_annot_blastn
###########################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

db = snakemake.input.db[:-4]

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_annot_blastn \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

text = "Filtering taxonomy tree on node "+str(snakemake.params.taxid)+" (named by user as "+snakemake.wildcards.tax+")"
f = open(snakemake.log.run, 'at')
f.write("## INFO: "+text+"\n")
f.close()

command = "$(which time) get_species_taxids.sh -t "+str(snakemake.params.taxid)+" > "+snakemake.params.tax_list+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) blastdbcmd -db "+db+" -taxidlist "+snakemake.params.tax_list+"|grep '>'|wc -l"
text = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## INFO: Resulting in "+text+" sequences from database\n")
f.close()

# Blast transcripts against the database to identify potential orthologs
command = "$(which time) blastn -task megablast"+\
          " -query " + snakemake.input.fa + \
          " -db " + db +\
          " -taxidlist "+snakemake.params.tax_list+\
          " -num_threads " + str(snakemake.threads) +\
          " -max_target_seqs 10"+\
          " -evalue 1e-5"+\
          " -mt_mode 1"+\
          " -outfmt '6 qaccver qlen qstart qend saccver slen sstart send length bitscore evalue pident qcovhsp mismatch gapopen staxid'"+\
          " -out " + snakemake.params.blast_full + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# Default order from BLAST outfmt6:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

command = '$(which time) Rscript -e  \'data.table::fwrite(data.table::fread("'+snakemake.params.blast_full+'")'+\
          '[V11<='+str(snakemake.params.eval_cutof)+' & V12>='+str(snakemake.params.pident_cutof)+\
          ' & (V9-V14)/V2>='+str(snakemake.params.cov_cutof)+',.SD[order(-V10,V11,V14)][1],by=V1][,.(V1,V5,V12,V9,V14,V15,V3,V4,V7,V8,V11,V10)]'+\
          ',"'+snakemake.output.blast_out+'",sep="\\\\t",quote=F,col.names=F,row.names=F)\' >> '+snakemake.log.run+' 2>&1'
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
