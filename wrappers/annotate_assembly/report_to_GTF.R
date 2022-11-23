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

map_to_uniprot_name = function(id, srcs) {
  print(paste0("looking for ",id," in ",paste0(srcs,collapse = ", ")))
  # id = "CR388227"
  # srcs = c("EMBL")
  res = data.table()
  for(src in srcs) {
    res = cbind(
      res,
      suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,src,".gz | grep -wF '",id,"'"), header = F, col.names = c("uniprot_acc","query"), key = "uniprot_acc")))
  }
  if(res[,.N]>0) {
    names = uniprot_gene_name[uniprot_acc %in% res$uniprot_acc]
    # names = suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"Gene_Name.gz | grep -wF '",res$uniprot_acc,"'"), header = F, col.names = c("uniprot_acc","uniprot_name"), key="uniprot_acc"))
    if(names[,.N]>0) {
      res = names[res, on="uniprot_acc"]
    } else {
      # not matching UniProt name found for UniProt accession
      return("")
    }
    if(res[,.N]>1) {
      # res = suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"UniProtKB-ID.gz | grep -wFf ",gfile), header = F, col.names = c("uniprot_acc","uniprot_id"), key="uniprot_acc"))[res, on="uniprot_acc"]
      # res = suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"HGNC.gz | grep -wFf ",gfile), header = F, col.names = c("uniprot_acc","hgnc_id"), key="uniprot_acc"))[res, on="uniprot_acc"]
      res[,uniprot_rev:=(function(acc) {
        ttab = suppressWarnings(fread(paste0('https://www.uniprot.org/uniprot/',uniprot_acc,'.txt'), sep="\n", header = F))
        return(ifelse(ttab[,.N]>0, ttab[V1 %like% "^ID",!grepl('Unreviewed',V1,fixed = T)], NA))
      })(uniprot_acc),by=uniprot_acc]
      ### Unsuccessful attempt to group hits by CRC hash, unfortunately, good hits are individual/non-redundant and we aim for those
      # res = res[suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"CRC64.gz | grep -wFf ",gfile), header = F, col.names = c("uniprot_acc","uniprot_hash"), key="uniprot_acc")), on="uniprot_acc"]
      # res[,.(uniprot_id_name=sub('_.*$','',uniprot_id), uniprot_id_org=sub('^[^_]_','',uniprot_id))][,.N,by=uniprot_id_name][order(N)]
      if(res[uniprot_rev == T & !is.na(uniprot_name), .N]>0) {
        res = res[uniprot_rev == T & !is.na(uniprot_name)]
      }
    }
    if(res[!is.na(uniprot_name), .N]>0) {
      return(res[!is.na(uniprot_name), .N, by=uniprot_name][which.max(N),uniprot_name])
    }
  }
  return("")
}

process_blast_evidence = function(dt, col) {
  # Take only the evidence with the lowest E-value
  best_ev = dt[get(col) != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(get(col),'`')),'^', simplify = T))[,
                                              V5:=as.numeric(sub('^E:','',V5))][which.min(V5), 
                                                      c(V1,paste0(V1,',',V3,',',V4,',E:',V5))]]
  if(grepl("refseq", col)) {
    # if RefSeq DB result
    source = "RefSeq"
    dbs = c("RefSeq", "RefSeq_NT")
  } else if(grepl("merops", col)) {
    # if MEROPS DB result
    source = "MEROPS"
    dbs = c("MEROPS")
    best_ev[1] = merops[id == best_ev[1], acc]
  } else if(grepl("_nt_",col)) {
    # if nucleotide DB result, take the first part of column name and paste it at the end of NCBI_nucleotide_XXX source string
    source = paste0("NCBI_Nucleotide_",sub('^(\\w)','\\U\\1',strsplit(col, '_')[[1]][1],perl = T))
    if(grepl('_',best_ev[1])) {
      dbs = c("RefSeq_NT")
    } else {
      best_ev[1] = sub('\\.[0-9]+$','',best_ev[1])
      dbs = c("EMBL")
    }
  } else {
    # if anything else just take the first part of column name as a source
    source = strsplit(col, '_')[[1]][1]
    dbs = c()
  }
  if(length(dbs)>0) {
    uname = map_to_uniprot_name(best_ev[1], dbs)
    if(uname != '') {
      usrc = "UniProt"
      return(c(uname, usrc, best_ev[2]))
    }
  }
  return(c(best_ev[1],source,best_ev[2]))
}

best_evidence = function(dt) {
  if(sum(dt$Kegg != '.')>0) {
    # split all KEGG evidence, take the first one, try to cross-reference against UniProt and return the whole line as gene_evidence
    best_ev = dt[Kegg != '.', .SD[1]$Kegg]
    uname = map_to_uniprot_name(sub("^KEGG:","",best_ev), "KEGG")
    if(uname != "") {
      usrc = "UniProt"
    } else {
      uname = sub("^KEGG:","",best_ev)
      usrc = "KEGG"
    }
    return(c(uname, usrc, best_ev))
  }
  if(sum(dt$sprot_Top_BLASTX_hit != '.')>0) {
    # split all sprot_Top_BLASTX_hit evidence, take the one with lowest E-value (V5), cross-reference against UniProt gene names and return the whole line as gene_evidence
    best_ev = dt[sprot_Top_BLASTX_hit != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(sprot_Top_BLASTX_hit,'`')),'^', simplify = T))[,
                                                            V5:=as.numeric(sub('^E:','',V5))][which.min(V5), 
                                                                                              c(V1,"UniProtKB/Swiss-Prot",paste0(V1,',',V3,',',V4,',E:',V5))]]
    uname = map_to_uniprot_name(best_ev[1], "UniProtKB-ID")
    if(uname != '') {
      usrc = "UniProt"
    } else {
      uname = sub('_[^_]+$','',best_ev[1])
      usrc = best_ev[2]
    }
    return(c(uname, usrc, best_ev[3]))
  }
  if(sum(dt$sprot_Top_BLASTP_hit != '.')>0) {
    # split all sprot_Top_BLASTP_hit evidence, take the one with lowest E-value (V5) and return name as gene_name as well as the whole line as gene_evidence
    best_ev = dt[sprot_Top_BLASTP_hit != '.', as.data.table(stri_split_fixed(unlist(stri_split_fixed(sprot_Top_BLASTP_hit,'`')),'^', simplify = T))[,
                                                            V5:=as.numeric(sub('^E:','',V5))][which.min(V5), 
                                                                                              c(V1,"UniProtKB/Swiss-Prot",paste0(V1,',',V3,',',V4,',E:',V5))]]
    uname = map_to_uniprot_name(best_ev[1], "UniProtKB-ID")
    if(uname != '') {
      usrc = "UniProt"
    } else {
      uname = sub('_[^_]+$','',best_ev[1])
      usrc = best_ev[2]
    }
    return(c(uname, usrc, best_ev[3]))
  }
  for(i in custom_prot_colnames) {
    if(sum(dt[[i]] != '.')>0) {
      # process all custom BLAST evidence one by one, split the first found and take the one with lowest E-value (V5), map it against UniProt and return ID as gene_name as well as the whole line as gene_evidence
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
  for(i in custom_nt_colnames) {
    if(sum(dt[[i]] != '.')>0) {
      # process all custom BLAST evidence one by one, split the first found and take the one with lowest E-value (V5), map it against UniProt and return ID as gene_name as well as the whole line as gene_evidence
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
  if(sum(dt$Pfam != '.')>0) {
    # split all Pfam evidence, take the one with lowest E-value (V5) and return domain name as gene_name as well as the whole line as gene_evidence
    return(dt[Pfam != '.', 
              as.data.table(stri_split_fixed(unlist(stri_split_fixed(Pfam,'`')),'^', simplify = T))[,
                            V5:=as.numeric(sub('^E:','',V5))][which.min(V5), 
                                    c(V2,"Pfam",paste0(.SD,collapse = ","))]])
  }
  if(sum(dt$SignalP != '.')>0) {
    # split all SignalP evidence, take the one with highest score (V3) and return unknown_protein as gene_name as well as the whole line as gene_evidence
    return(dt[SignalP != '.', 
              as.data.table(stri_split_fixed(unlist(stri_split_fixed(SignalP,'`')),'^', simplify = T))[which.max(V3), 
                             c("unknown_protein","SignalP",paste0(V1,'-',V2,',',V3))]])
  }
}

# compress_evidence = function(vec) {
#   if(any(vec!='.')) {
#     return(paste0(vec[which(vec!='.')],collapse = '`'))
#   } else {
#     return('.')
#   }
# }

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

path_to_uniprot_mapping = "/mnt/ssd/ssd_1/workspace/martin/uniprot_id_mapping/"
# setwd("/mnt/ssd/ssd_1/snakemake/")
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
merops = fread("/mnt/ssd/ssd_3/references/general/MEROPS_DB/merops_pepunit.info.tsv", sep = '\t', header = F, col.names = c("id","desc", "org", "acc", "sub", "type", "pos", "name"))
print("loading uniprot gene_name reference")
uniprot_gene_name = fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"Gene_Name.gz"), header = F, sep = '\t', col.names = c("uniprot_acc","uniprot_name"), key = 'uniprot_acc')
# embl_to_uniprot = suppressWarnings(fread(cmd=paste0("zcat ",path_to_uniprot_mapping,"EMBL.gz"), header = F, sep = '\t', col.names = c("uniprot_acc","query"), key = 'query'))

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
  custom_prot_colnames = c(custom_prot_colnames, tab[,colnames(.SD),.SDcols=patterns(i)])
}
custom_nt_colnames = c()
for(i in custom_nt_taxid) {
  custom_nt_colnames = c(custom_nt_colnames, tab[,colnames(.SD),.SDcols=patterns(i)])
}

gtf = tab[prot_coords != "."][,{
  print(paste0(.BY,' (',.I,'/',.NGRP,')'))
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
    gene_name = "gene_with_evidence" #hit[1]
    gene_source = "see_transcript" #hit[2]
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
  }, by=transcript_id]

print("writing final annotation")
fwrite(gtf[,.SD,.SDcols=!c("prot_id")],
       output_name,
       sep = '\t',
       col.names = F,
       row.names = F,
       quote = F)
