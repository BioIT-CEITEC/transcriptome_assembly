library(data.table)
library(stringi)

sum_evidence = function(dt) {
  dt[, sum(
    # sum(Kegg != '.'),
    sum(sprot_Top_BLASTX_hit != '.'),
    sum(sprot_Top_BLASTP_hit != '.'),
    sapply(custom_prot_colnames, function(x) sum(get(x) != '.')),
    sapply(custom_nt_colnames, function(x) sum(get(x) != '.')),
    sum(Pfam != '.'),
    sum(SignalP != '.'),
    sum(RNAMMER != '.')
  )]
}

process_blast_evidence = function(dt, col) {
  # Take only the evidence with the lowest E-value
  best_ev = dt[get(col) != '.', 
               as.data.table(
                 stri_split_fixed(
                   unlist(stri_split_fixed(get(col),'`')),
                   '^', simplify = T))[,
                                       V5:=as.numeric(sub('^E:','',V5))][which.min(V5),
                                                                         c(V1,paste0(V1,',',V3,',',V4,',E:',V5))]]
  if(grepl("refseq", col)) {
    # if RefSeq DB result
    source = "RefSeq"
  } else if(grepl("merops", col)) {
    # if MEROPS DB result
    source = "MEROPS"
    best_ev[1] = merops[id == best_ev[1], acc]
  } else if(grepl("_nt_",col)) {
    # if nucleotide DB result, take the first part of column name and paste it at the end of NCBI_nucleotide_XXX source string
    source = paste0("NCBI_Nucleotide_",sub('^(\\w)','\\U\\1',strsplit(col, '_')[[1]][1],perl = T))
    if(grepl('_',best_ev[1])) {
      source = paste0(source,"(RefSeq_id)")
    } else {
      best_ev[1] = sub('\\.[0-9]+$','',best_ev[1])
      source = paste0(source,"(EMBL_id)")
    }
  } else {
    # if anything else just take the first part of column name as a source
    source = strsplit(col, '_')[[1]][1]
  }
  return(c(best_ev[1],source,best_ev[2]))
}

best_evidence = function(dt) {
  if(sum(dt$Kegg != '.')>0) {
    # split all KEGG evidence, take the first one, try to cross-reference against UniProt and return the whole line as gene_evidence
    best_ev = dt[Kegg != '.', .SD[1]$Kegg]
    uname = sub("^KEGG:","",best_ev)
    usrc = "KEGG"
    return(c(uname, usrc, best_ev))
  }
  if(sum(dt$sprot_Top_BLASTX_hit != '.')>0) {
    # split all sprot_Top_BLASTX_hit evidence, take the one with lowest E-value (V5) and return the whole line as gene_evidence
    best_ev = dt[sprot_Top_BLASTX_hit != '.', 
                 as.data.table(
                   stri_split_fixed(
                     unlist(stri_split_fixed(sprot_Top_BLASTX_hit,'`')),
                     '^', simplify = T))[,
                                         V5:=as.numeric(sub('^E:','',V5))][which.min(V5),
                                                                           c(stri_split_fixed(V1,"|")[[1]][3],
                                                                             "UniProtKB/Swiss-Prot",
                                                                             paste0(V1,',',V3,',',V4,',E:',V5))]]
    return(best_ev)
  }
  if(sum(dt$sprot_Top_BLASTP_hit != '.')>0) {
    # split all sprot_Top_BLASTP_hit evidence, take the one with lowest E-value (V5) and return name as gene_name as well as the whole line as gene_evidence
    best_ev = dt[sprot_Top_BLASTP_hit != '.', 
                 as.data.table(
                   stri_split_fixed(
                     unlist(stri_split_fixed(sprot_Top_BLASTP_hit,'`'))
                     ,'^', simplify = T))[,
                                          V5:=as.numeric(sub('^E:','',V5))][which.min(V5), 
                                                                            c(stri_split_fixed(V1,"|")[[1]][3],
                                                                              "UniProtKB/Swiss-Prot",
                                                                              paste0(V1,',',V3,',',V4,',E:',V5))]]
    return(best_ev)
  }
  for(i in custom_prot_colnames) {
    if(sum(dt[[i]] != '.')>0) {
      # process all custom BLAST evidence one by one, split the first found and take the one with lowest E-value (V5), map it against UniProt and return ID as gene_name as well as the whole line as gene_evidence
      return(process_blast_evidence(dt, i))
    }
  }
  for(i in custom_nt_colnames) {
    if(sum(dt[[i]] != '.')>0) {
      # process all custom BLAST evidence one by one, split the first found and take the one with lowest E-value (V5), map it against UniProt and return ID as gene_name as well as the whole line as gene_evidence
      return(process_blast_evidence(dt, i))
    }
  }
  if(sum(dt$Pfam != '.')>0) {
    # split all Pfam evidence, take the one with lowest E-value (V5) and return domain name as gene_name as well as the whole line as gene_evidence
    return(dt[Pfam != '.', 
              as.data.table(
                stri_split_fixed(
                  unlist(stri_split_fixed(Pfam,'`'))
                  ,'^', simplify = T))[,
                                       V5:=as.numeric(sub('^E:','',V5))][which.min(V5),
                                                                         c(V2,"Pfam",paste0(.SD,collapse = ","))]])
  }
  if(sum(dt$SignalP != '.')>0) {
    # split all SignalP evidence, take the one with highest score (V3) and return unknown_protein as gene_name as well as the whole line as gene_evidence
    return(dt[SignalP != '.', 
              as.data.table(
                stri_split_fixed(
                  unlist(stri_split_fixed(SignalP,'`')),
                  '^', simplify = T))[which.max(V3),
                                      c("unknown_protein","SignalP",paste0(V1,'-',V2,',',V3))]])
  }
}

path_to_uniprot_mapping = "/mnt/ssd/ssd_1/workspace/martin/uniprot_id_mapping/"
# setwd("/mnt/ssd/ssd_1/sequia/6744__transcriptome_assembly__genome_guided_TA__221025/")
args = commandArgs(trailingOnly=TRUE)
# args = c("results/annotation/trinotate_report_wt_seqs.xls",
#          "results/annotation/transdecoder_dir/merged_samples.merged.fa.transdecoder.pep",
#          "results/annotation/merged_samples.trinotate_annotation.gtf",
#          "merops",
#          "arabidopsis,brassicaceae")
input_name  = args[1]
trans_input = args[2]
output_name = args[3]
custom_prot_db = strsplit(args[4],',')[[1]]
custom_nt_taxid= strsplit(args[5],',')[[1]]

print("loading merops info file")
merops = fread("/mnt/shared/CFBioinformatics/references_backup/general/MEROPS_DB/merops_pepunit.info.tsv", sep = '\t', header = F, col.names = c("id","desc", "org", "acc", "sub", "type", "pos", "name"))
# print("loading uniprot gene_name reference")
# uniprot_gene_name = fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"Gene_Name.gz"), header = F, sep = '\t', col.names = c("uniprot_acc","uniprot_name"), key = 'uniprot_acc')
# embl_to_uniprot = suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"EMBL.gz"), header = F, sep = '\t', col.names = c("uniprot_acc","query"), key = 'query'))
# refseq_nt_to_uniprot = suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"RefSeq_NT.gz"), header = F, sep = '\t', col.names = c("uniprot_acc","query"), key = 'query'))

print("loading transdecoder.pep")
trans = fread(cmd=paste0('grep "^>" ',trans_input), sep=" ", header = F)
trans[,V1:=sub('^>','',V1)]
if(trans[!V4 %like% "type:",.N]>0 || trans[!V6 %like% "score=",.N]>0) {
  cat("Warning: Some header lines from ",trans_input," have incorrect format! Either ORF type or score is missing or displaced!",
      trans[1:3, apply(.SD, 1, paste0, collapse = " ")], sep = "\n")
}
trans[,score:=as.numeric(sub('score=','',tstrsplit(V6,',')[[2]]))]
trans[,ORF_type:=sub('type:','',V4)]

print("loading trinity report")
tab = fread(input_name)
tab[,N:=.N,by=transcript_id]
tab[,trans_len:=nchar(transcript)]
# colnames(tab)
# [1] "#gene_id"                 "transcript_id"            "sprot_Top_BLASTX_hit"     "RNAMMER"                  "prot_id"                 
# [6] "prot_coords"              "sprot_Top_BLASTP_hit"     "merops_BLASTX"            "refseq_BLASTX"            "actinopterygii_nt_BLASTX"
# [11] "cypriniformes_nt_BLASTX"  "merops_BLASTP"            "refseq_BLASTP"            "Pfam"                     "SignalP"                 
# [16] "TmHMM"                    "eggnog"                   "Kegg"                     "gene_ontology_BLASTX"     "gene_ontology_BLASTP"    
# [21] "gene_ontology_Pfam"       "transcript"               "peptide"                  "N"                        "trans_len" 

tab[prot_coords != ".", start  := as.numeric(gsub('^([0-9]+).*', '\\1', prot_coords, perl=T))]
tab[prot_coords != ".", end    := as.numeric(gsub('^[0-9]+-([0-9]+).*', '\\1', prot_coords, perl=T))]
tab[prot_coords != ".", strand := gsub('.*\\[([+-])\\]', '\\1', prot_coords, perl=T)]

custom_prot_colnames = c()
for(i in custom_prot_db) {
  #custom_prot_colnames = c(custom_prot_colnames, tab[,colnames(.SD),.SDcols=patterns(i)])
  new_names = grep(i, colnames(tab), value=T)
  if(length(new_names)==0) {
    print(paste("ID:",i,"wasn't found in colnames of input:",input_name))
  } else {
    custom_prot_colnames = c(custom_prot_colnames, new_names)
  }
}
custom_nt_colnames = c()
for(i in custom_nt_taxid) {
  #custom_nt_colnames = c(custom_nt_colnames, tab[,colnames(.SD),.SDcols=patterns(i)])
  new_names = grep(i, colnames(tab), value=T)
  if(length(new_names)==0) {
    print(paste("ID:",i,"wasn't found in colnames of input:",input_name))
  } else {
    custom_nt_colnames = c(custom_nt_colnames, new_names)
  }
}

setnames(tab, '#gene_id', "gene_id")
print("extracting annotation info")
start_time <- Sys.time()
annot_tab = tab[prot_coords != ".", {
  evid = sum_evidence(.SD)
  if(evid == 0) {
    gene_name = "unknown_gene"
    gene_source = "TransDecoder"
    gene_evid = "ORF_prediction"
    gene_info = "0_evidence_sources"
  } else {
    hit = best_evidence(.SD)
    gene_name = hit[1]
    gene_source = hit[2]
    gene_evid = hit[3]
    gene_info = ifelse(evid==1, paste0(evid,"_evidence_source"), paste0("best_evidence_of_",evid,"_sources"))
  }
  .SD[order(start),
    list(
      chr_id=gene_id,
      transcript_id=.BY,
      gene_id=.BY,
      gene_start=start,
      gene_end=end,
      gene_score=trans[V1==.BY, score],
      gene_strand=strand,
      gene_phase="0",
      gene_name=gene_name,
      gene_biotype="protein_coding",
      gene_source=gene_source,
      gene_evidence=gene_evid,
      gene_orf=trans[V1==.BY, ORF_type],
      gene_info=gene_info)
  ]
}, by=prot_id]
print(Sys.time() - start_time)

print("formatting final annotation")
start_time <- Sys.time()
gtf = annot_tab[order(chr_id, gene_start),{
  # Line with gene information
  gene_attrs = paste0('gene_id "',gene_id,'";',
                      ' gene_name "',gene_name[1],'";',
                      ' gene_biotype "',gene_biotype[1],'";',
                      ' gene_source "',gene_source[1],'";',
                      ' gene_evidence "',gene_evidence[1],'";',
                      ' gene_info "',gene_info[1],'";')
  # Line with transcript information
  trans_attrs = paste0(gene_attrs,
                       ' transcript_id "',transcript_id,'";',
                       ' transcript_orf "',gene_orf[1],'";')
  # Line with exon information
  exon_attrs = paste0(trans_attrs,
                      ' exon_id "exon.',transcript_id,'";',
                      ' exon_number "1";')
  gene_line = list(
    V1=rep(chr_id, 3),
    V2=rep("Trinotate", 3),
    V3=c("gene", "transcript", "exon"),
    V4=rep(gene_start, 3),
    V5=rep(gene_end, 3),
    V6=rep(gene_score, 3),
    V7=rep(gene_strand, 3),
    V8=rep(gene_phase, 3),
    V9=c(gene_attrs, trans_attrs, exon_attrs))
}, by=.(chr_id,prot_id)][,.SD,.SDcols=!c("prot_id","chr_id")]
print(Sys.time() - start_time)

print("writing final annotation")
fwrite(gtf,
       output_name,
       sep = '\t',
       col.names = F,
       row.names = F,
       quote = F)
