###########################################
# wrapper for rule: annotate_assembly
###########################################
import os
import sys
import math
import subprocess
import uuid
from snakemake.shell import shell

shell.executable("/bin/bash")

trinotate = "Trinotate"

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: annotate_assembly \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

trin_report = snakemake.output.trin_rep
trin_report_full = snakemake.params.trin_rep_full
# trin_report_eval = snakemake.output.trin_rep_eval
# trin_GO_annot = snakemake.output.go_rep
#sqlite = snakemake.params.sqlite
sqlite = os.path.join(snakemake.params.tmp, "trinotate_sqlite_"+uuid.uuid4().hex, os.path.basename(snakemake.params.sqlite))

command = "mkdir -p "+os.path.dirname(sqlite)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cp "+snakemake.input.sqlite+" "+sqlite+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# Build the boilerplate trinotate sqlite database
command = "$(which time) "+trinotate + " " + sqlite + " init --gene_trans_map " + snakemake.input.gtmap + " " + \
          " --transcript_fasta " + snakemake.input.fa + " --transdecoder_pep " + snakemake.input.orfs + " " + \
          " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# Load annotation results into the SQLite database
command = "$(which time) "+trinotate + " " + sqlite + " LOAD_swissprot_blastp " + snakemake.input.blastp_uniprot + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) "+trinotate + " " + sqlite + " LOAD_swissprot_blastx " + snakemake.input.blastx_uniprot + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) "+trinotate + " " + sqlite + " LOAD_pfam " + snakemake.input.pfam + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) "+trinotate + " " + sqlite + " LOAD_signalp " + snakemake.input.signalp + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

for inp in snakemake.input.custom:
  ins = os.path.basename(inp).split("_")
  prog = ins[0] if ins[0] != "blastn" else "blastp"
  name = ins[1] if ins[1] != "nt" else ins[3].split(".")[0]+"_nt"
  print("Loading custom: "+name+" as prog "+prog)

  # --outfmt6 <file> --prog <blastp|blastx> --dbtype <database_name>
  command = "$(which time) "+trinotate + " " + sqlite + " LOAD_custom_blast"+\
            " --outfmt6 " + inp +\
            " --prog "+prog+\
            " --dbtype "+name+\
            " >> " + snakemake.log.run + " 2>&1 "
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

# MAYBE for the future
# Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
# Trinotate Trinotate.sqlite LOAD_rnammer $INPUT_FILE.rnammer.gff

# Use trinotate to import transcript name information into the report
command = "$(which time) import_transcript_names.pl " + sqlite + " " + snakemake.input.gtmap + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# Generate annotation reports
command = "$(which time) "+trinotate + " " + sqlite + " report > " + trin_report + " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) "+trinotate + " " + sqlite + " report --incl_pep --incl_trans > " + trin_report_full + " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# # Trinotate Trinotate.sqlite report -E 0.00001 > $INPUT_FILE.trinotate_annotation_report_evalue00001.xls
# command = "$(which time) "+trinotate + " " + sqlite + " report -E 0.00001 > " + trin_report_eval + " 2>> " + snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# # Extract GO terms from the report
# # parameter -I/--include_ancestral_terms is not working because of missing obo/go-basic.obo.gz
# #command = "extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls " + trin_report + " -G -I > " + snakemake.output.go_rep + " 2>> " + snakemake.log.run
# command = "$(which time) extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls " + trin_report + " -G > " + snakemake.output.go_rep + " 2>> " + snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

command = "cp "+sqlite+" "+snakemake.params.sqlite+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) Rscript "+snakemake.params.rscript+\
          " "+trin_report_full+\
          " "+snakemake.input.orfs+\
          " "+snakemake.output.gtf+\
          " "+",".join(snakemake.params.custom_prot_db)+\
          " "+",".join(snakemake.params.custom_nt_taxid)+\
          " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
