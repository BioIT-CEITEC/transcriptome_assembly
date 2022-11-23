library(data.table)
library(stringi)
setwd("/mnt/ssd/ssd_1/snakemake/")

sum_evidence = function(dt) {
  dt[, sum(
          # sum(Kegg != '.'),
          sum(sprot_Top_BLASTX_hit != '.'),
          sum(sprot_Top_BLASTP_hit != '.'),
          sapply(custom_prot_colnames, function(x) sum(get(x) != '.')),
          sum(Pfam != '.'),
          sapply(custom_nt_colnames, function(x) sum(get(x) != '.')),
          sum(SignalP != '.'),
          sum(RNAMMER != '.')
          )]
}

process_blast_evidence = function(dt, col) {
  if(grepl("refseq", col)) {
    # if RefSeq DB result
    source = "RefSeq"
  } else if(grepl("merops", col)) {
    # if MEROPS DB result
    source = "MEROPS"
  } else if(grepl("_nt_",col)) {
    # if nucleotide DB result, take the first part of column name and paste it at the end of NCBI_nucleotide_XXX source string
    source = paste0("NCBI_Nucleotide_",sub('^(\\w)','\\U\\1',strsplit(col, '_')[[1]][1],perl = T))
  } else {
    # if anything else just take the first part of column name as a source
    source = strsplit(col, '_')[[1]][1]
  }
  dt[get(col) != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(get(col),'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,source,paste0(V1,',',V3,',',V4,',E:',V5))]]
}

best_evidence = function(dt) {
  if(sum(dt$Kegg != '.')>0) {
    # split all KEGG evidence, take the first one and return whole line as gene_name as well as gene_evidence
    return(dt[Kegg != '.', c(.SD[1]$Kegg,"KEGG",.SD[1]$Kegg)])
  }
  if(sum(dt$sprot_Top_BLASTX_hit != '.')>0) {
    # split all sprot_Top_BLASTX_hit evidence, take the one with lowest E-value (V5) and return name as gene_name as well as the whole line as gene_evidence
    return(dt[sprot_Top_BLASTX_hit != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(sprot_Top_BLASTX_hit,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(sub('_[^_]+$','',V1),"UniProtKB/Swiss-Prot",paste0(V1,',',V3,',',V4,',E:',V5))]])
  }
  if(sum(dt$sprot_Top_BLASTP_hit != '.')>0) {
    # split all sprot_Top_BLASTP_hit evidence, take the one with lowest E-value (V5) and return name as gene_name as well as the whole line as gene_evidence
    return(dt[sprot_Top_BLASTP_hit != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(sprot_Top_BLASTP_hit,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(sub('_[^_]+$','',V1),"UniProtKB/Swiss-Prot",paste0(V1,',',V3,',',V4,',E:',V5))]])
  }
  for(i in custom_prot_colnames) {
    if(sum(dt[[i]] != '.')>0) {
      # process all custom BLAST evidence one by one, split the first found and take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
      return(process_blast_evidence(dt, i))
    }
  }
  # if(sum(dt$refseq_BLASTX != '.')>0) {
  #   # split all refseq_BLASTX evidence, take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
  #   return(dt[refseq_BLASTX != '.',as.data.table(stri_split_fixed(unlist(stri_split_fixed(refseq_BLASTX,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,"RefSeq",paste0(V1,',',V3,',',V4,',E:',V5))]])
  # }
  # if(sum(dt$refseq_BLASTP != '.')>0) {
  #   # split all refseq_BLASTP evidence, take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
  #   return(dt[refseq_BLASTP != '.',as.data.table(stri_split_fixed(unlist(stri_split_fixed(refseq_BLASTP,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,"RefSeq",paste0(V1,',',V3,',',V4,',E:',V5))]])
  # }
  # if(sum(dt$merops_BLASTX != '.')>0) {
  #   # split all merops_BLASTX evidence, take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
  #   return(dt[merops_BLASTX != '.',as.data.table(stri_split_fixed(unlist(stri_split_fixed(merops_BLASTX,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,"MEROPS",paste0(V1,',',V3,',',V4,',E:',V5))]])
  # }
  # if(sum(dt$merops_BLASTP != '.')>0) {
  #   # split all merops_BLASTP evidence, take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
  #   return(dt[merops_BLASTP != '.',as.data.table(stri_split_fixed(unlist(stri_split_fixed(merops_BLASTP,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,"MEROPS",paste0(V1,',',V3,',',V4,',E:',V5))]])
  # }
  if(sum(dt$Pfam != '.')>0) {
    # split all Pfam evidence, take the one with lowest E-value (V5) and return domain name as gene_name as well as the whole line as gene_evidence
    return(dt[Pfam != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(Pfam,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V2,"Pfam",paste0(.SD,collapse = ","))]])
  }
  for(i in custom_nt_colnames) {
    if(sum(dt[[i]] != '.')>0) {
      # process all custom BLAST evidence one by one, split the first found and take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
      return(process_blast_evidence(dt, i))
    }
  }
  # if(sum(dt$actinopterygii_nt_BLASTX != '.')>0) {
  #   # split all actinopterygii_nt_BLASTX evidence, take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
  #   return(dt[actinopterygii_nt_BLASTX != '.',as.data.table(stri_split_fixed(unlist(stri_split_fixed(actinopterygii_nt_BLASTX,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,"NCBI_Nucleotide_Actinopterygii",paste0(V1,',',V3,',',V4,',E:',V5))]])
  # }
  # if(sum(dt$cypriniformes_nt_BLASTX != '.')>0) {
  #   # split all cypriniformes_nt_BLASTX evidence, take the one with lowest E-value (V5) and return ID as gene_name as well as the whole line as gene_evidence
  #   return(dt[cypriniformes_nt_BLASTX != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(cypriniformes_nt_BLASTX,'`')),'^', simplify = T))[,V5:=as.numeric(sub('^E:','',V5))][which.min(V5), c(V1,"NCBI_Nucleotide_Cypriniformes",paste0(V1,',',V3,',',V4,',E:',V5))]])
  # }
  if(sum(dt$SignalP != '.')>0) {
    # split all SignalP evidence, take the one with highest score (V3) and return unknown_protein as gene_name as well as the whole line as gene_evidence
    return(dt[SignalP != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(SignalP,'`')),'^', simplify = T))[which.max(V3), c("unknown_protein","SignalP",paste0(V1,'-',V2,',',V3))]])
  }
}

compress_evidence = function(vec) {
  if(any(vec!='.')) {
    return(paste0(vec[which(vec!='.')],collapse = '`'))
  } else {
    return('.')
  }
}

make_transcript_and_exon = function(dt, attrs) {
  dt[order(start),{
    # Line with transcript information
    trans_attrs = paste0(attrs,
                   ' transcript_id "',.BY,'";',
                   ' transcript_orf "',trans[V1==.BY, ORF_type],'";')
    evid = sum_evidence(.SD)
    if(evid == 0) {
      trans_attrs = paste0(trans_attrs,
                           ' transcript_info "0_evidence_sources";')
    } else {
      hit = best_evidence(.SD)
      trans_name = hit[1]
      trans_source = hit[2]
      trans_evid = hit[3]
      trans_info = ifelse(evid==1, paste0(evid,"_evidence_source"), paste0("best_evidence_of_",evid,"_sources"))
      trans_attrs = paste0(trans_attrs,
                     ' transcript_name "',trans_name,'";',
                     ' transcript_source "',trans_source,'";',
                     ' transcript_evidence "',trans_evid,'";',
                     ' transcript_info "',trans_info,'";')
    }
    # Line with exon information
    exon_attrs = paste0(trans_attrs,
                        ' exon_id "exon.',.BY,'";',
                        ' exon_number "1";')
    list(
      V2=rep("Trinotate",2), 
      V3=c("transcript","exon"), 
      V4=rep(start,2), 
      V5=rep(end,2), 
      V6=rep(trans[V1==.BY, score],2),
      V7=rep(strand,2), 
      V8=rep("0",2), 
      V9=c(trans_attrs,exon_attrs)) 
  }, prot_id]
}


abramis = T
if(abramis) {
  # Abramis
  trans_input = "stage520_Abramis_rutilus_TAssembly.Abramis_brama_cejn/denovo_transcriptome_assembly/assembly/abramis/annotation/transdecoder_dir/abramis.merged.fa.transdecoder.pep"
  input_name  = "stage520_Abramis_rutilus_TAssembly.Abramis_brama_cejn/denovo_transcriptome_assembly/assembly/abramis/annotation/trinotate_report_wt_seqs.xls"
} else {
  # Rutilus
  trans_input = "stage519_Abramis_rutilus_TAssembly.Rutilus_rutilus_plotice/denovo_transcriptome_assembly/assembly/rutilus/annotation/transdecoder_dir/rutilus.merged.fa.transdecoder.pep"
  input_name  = "stage519_Abramis_rutilus_TAssembly.Rutilus_rutilus_plotice/denovo_transcriptome_assembly/assembly/rutilus/annotation/trinotate_report_wt_seqs.xls"
}
custom_prot_db = c("refseq","merops")
custom_nt_taxid= c("actinopterygii","cypriniformes")

trans = fread(cmd=paste0('grep "^>" ',trans_input), sep=" ", header = F)
trans[,V1:=sub('^>','',V1)]
if(trans[!V4 %like% "type:",.N]>0 || trans[!V6 %like% "score=",.N]>0) {
  cat("Warning: Some header lines from ",trans_input," have incorrect format! Either ORF type or score is missing or displaced!",
               trans[1:3, apply(.SD, 1, paste0, collapse = " ")], sep = "\n")
}
trans[,score:=as.numeric(sub('score=','',tstrsplit(V6,',')[[2]]))]
trans[,ORF_type:=sub('type:','',V4)]

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
  custom_prot_colnames = c(custom_prot_colnames, tab[,colnames(.SD),.SDcols=patterns(i)])
}
custom_nt_colnames = c()
for(i in custom_nt_taxid) {
  custom_nt_colnames = c(custom_nt_colnames, tab[,colnames(.SD),.SDcols=patterns(i)])
}

# gtf = tab[prot_coords != "." & transcript_id == "trinity_21.TRINITY_DN2347_c1_g2_i1", {
gtf = tab[prot_coords != ".", {
  # Line with gene information
  gene_strand = .SD[1]$strand
  evid = sum_evidence(.SD)
  if(evid == 0) {
    gene_name = paste0("unknown_gene")
    gene_source = "TransDecoder"
    gene_biotype = "protein_coding"
    gene_info = "ORF_prediction"
  } else {
    hit = best_evidence(.SD)
    gene_name = hit[1]
    gene_source = hit[2]
    gene_biotype = "protein_coding"
    gene_info = fifelse(evid==1,paste0(evid,"_evidence_source"),paste0("best_evidence_of_",evid,"_sources"))
  }
  attrs = paste0('gene_id "',.BY,'";',
                 ' gene_name "',gene_name,'";',
                 ' gene_source "',gene_source,'";',
                 ' gene_biotype "',gene_biotype,'";',
                 ' gene_info "',gene_info,'";')
  gene_line = list(
    V2="Trinotate", 
    V3="gene", 
    V4=1, 
    V5=unique(trans_len), 
    V6=".", 
    V7=gene_strand, 
    V8=".", 
    V9=attrs)
  
  rbind(gene_line,
        make_transcript_and_exon(.SD, attrs),
        fill=T)
  }, transcript_id]

fwrite(gtf[,.SD,.SDcols=!c("prot_id")],
       sub('xls$','gtf',input_name),
       sep = '\t',
       col.names = F,
       row.names = F,
       quote = F)

#####################################################
## This part wa used to merge overlapping exons
# tab = rbind(
#   tab[prot_coords!='.', {
#     .SD[order(start)][!(start <= data.table::shift(end,1,0)|end>=data.table::shift(start,1,end[.N],"lead")),
#                       .(transcript_id,
#                         sprot_Top_BLASTX_hit,
#                         RNAMMER,
#                         prot_coords,
#                         sprot_Top_BLASTP_hit,
#                         merops_BLASTX,
#                         refseq_BLASTX,
#                         actinopterygii_nt_BLASTX,
#                         cypriniformes_nt_BLASTX,
#                         merops_BLASTP,
#                         refseq_BLASTP,
#                         Pfam,
#                         SignalP,
#                         TmHMM,
#                         eggnog,
#                         Kegg,
#                         gene_ontology_BLASTX,
#                         gene_ontology_BLASTP,
#                         gene_ontology_Pfam,
#                         transcript,
#                         trans_len,
#                         start,
#                         end,
#                         strand)]
#   },transcript_id],
#   tab[prot_coords!='.', {
#     .SD[order(start)][start < data.table::shift(end,1,0)|end>data.table::shift(start,1,end[.N],"lead"),
#                       .(transcript_id=transcript_id[1],
#                         sprot_Top_BLASTX_hit=compress_evidence(sprot_Top_BLASTX_hit),
#                         RNAMMER=compress_evidence(RNAMMER),
#                         prot_coords=paste0(min(start),'-',max(end),'[',strand[1],']'),
#                         sprot_Top_BLASTP_hit=compress_evidence(sprot_Top_BLASTP_hit),
#                         merops_BLASTX=compress_evidence(merops_BLASTX),
#                         refseq_BLASTX=compress_evidence(refseq_BLASTX),
#                         actinopterygii_nt_BLASTX=compress_evidence(actinopterygii_nt_BLASTX),
#                         cypriniformes_nt_BLASTX=compress_evidence(cypriniformes_nt_BLASTX),
#                         merops_BLASTP=compress_evidence(merops_BLASTP),
#                         refseq_BLASTP=compress_evidence(refseq_BLASTP),
#                         Pfam=compress_evidence(Pfam),
#                         SignalP=compress_evidence(SignalP),
#                         TmHMM=compress_evidence(TmHMM),
#                         eggnog=compress_evidence(eggnog),
#                         Kegg=compress_evidence(Kegg),
#                         gene_ontology_BLASTX=compress_evidence(gene_ontology_BLASTX),
#                         gene_ontology_BLASTP=compress_evidence(gene_ontology_BLASTP),
#                         gene_ontology_Pfam=compress_evidence(gene_ontology_Pfam),
#                         transcript=transcript[1],
#                         trans_len=trans_len[1],
#                         start=min(start),
#                         end=max(end),
#                         strand=strand[1])]
#   },transcript_id])
# tab[prot_coords!='.', exon_overlap:=sum(.SD[order(start), start < data.table::shift(end,1,0)]),transcript_id]
# View(tab[exon_overlap>0])