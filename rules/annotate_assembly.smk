
def annotate_assembly_inputs(wc):
    inputs = {}
    inputs['fa'] = "results/assembly/merged/merged_samples.merged.fa"
    inputs['orfs'] = "results/annotation/transdecoder_dir/merged_samples.merged.fa.transdecoder.pep"
    inputs['gtmap'] = "results/assembly/merged/merged_samples.merged.fa.gene_trans_map"
    inputs['blastx_uniprot'] = "results/annotation/blastx_uniprot_processed.tsv"
    inputs['blastp_uniprot'] = "results/annotation/blastp_uniprot_processed.tsv"
    inputs['pfam'] = "results/annotation/TrinotatePFAM.out"
    inputs['pred_gff'] = "results/annotation/transdecoder_dir/merged_samples.merged.fa.transdecoder.gff3"
    inputs['signalp']  = "results/annotation/signalp.gff3"
    inputs['sqlite'] = os.path.join(GLOBAL_TMPD_PATH,"Trinotate.sqlite")
    inputs['custom'] = expand("results/annotation/blast{custom}.tsv", custom=["x_"+i+"_processed" for i in custom_prot_dbs]+["p_"+i+"_processed" for i in custom_prot_dbs]+["n_nt_processed.taxids_"+i.split(':')[0].lower() for i in nt_taxids])
    # TBD: Uniprot-mapping DBS: EMBL, MEROPS, RefSeq, RefSeq_NT
    return inputs

rule annotate_assembly:
    input:  unpack(annotate_assembly_inputs),
    output: trin_rep      = "results/annotation/trinotate_report.xls",
            trin_rep_full = "results/annotation/trinotate_report_wt_seqs.xls",
            gtf = "results/annotation/merged_samples.trinotate_annotation.gtf",
            # go_rep = "results/annotation/go_annotations.txt",
    log:    run = "logs/annotate_assembly.log",
    threads: 1
    resources:  mem = 90
    params: prefix = "results/annotation",
            sqlite = "results/annotation/Trinotate.sqlite",
            trin_rep_full = "results/annotation/trinotate_report_wt_seqs.xls",
            tmp = GLOBAL_TMPD_PATH,
            custom_prot_db = custom_prot_dbs,
            custom_nt_taxid= [i.split(':')[0].lower() for i in nt_taxids],
            rscript = workflow.basedir + "/wrappers/annotate_assembly/report_to_GTF.new_approach.R",
    conda:  "../wrappers/annotate_assembly/env.yaml"
    script: "../wrappers/annotate_assembly/script.py"
    
    
rule prepare_annot_signalp:
    input:  pep = "results/annotation/transdecoder_dir/longest_orfs.pep",
            # pep = "results/annotation/transdecoder_dir/merged_samples.merged.fa.transdecoder.pep",
            # licence = os.path.join(GLOBAL_REF_PATH,"general/SignalP/signalp-4.1g.Linux.tar.gz"),
            licence = os.path.join(GLOBAL_REF_PATH,"general/SignalP/signalp-5.0b.Linux.tar.gz"),
    output: ref = "results/annotation/signalp.gff3",
            fa  = "results/annotation/signalp.fa",
            tab = "results/annotation/signalp.out",
    log:    run = "logs/prepare_annot_signalp.log",
    threads: 30
    resources:  mem = 10
    params: prefix = "results/annotation/signalp",
            organism = "euk", #TODO: [euk, arch, gram+, gram-]
            tmp = GLOBAL_TMPD_PATH+"/",
            tool = "signalp5",
    conda:  "../wrappers/prepare_annot_signalp/env.yaml"
    script: "../wrappers/prepare_annot_signalp/script.py"


rule prepare_annot_transd_pred:
    input:  fa = "results/assembly/merged/merged_samples.merged.fa",
            pfam_out = "results/annotation/TrinotatePFAM.out",
            blastp_out = "results/annotation/blastp_uniprot_processed.tsv",
            orfs_pep = "results/annotation/transdecoder_dir/longest_orfs.pep",
    output: pred_out = "results/annotation/transdecoder_dir/merged_samples.merged.fa.transdecoder.gff3",
            pred_pep = "results/annotation/transdecoder_dir/merged_samples.merged.fa.transdecoder.pep",
    log:    run = "logs/prepare_annot_transd_pred.log"
    threads: 1
    resources:  mem = 30
    params: prefix = "results/annotation",
            transdecoder_dir = "results/annotation/transdecoder_dir/"
    conda:  "../wrappers/prepare_annot_transd_pred/env.yaml"
    script: "../wrappers/prepare_annot_transd_pred/script.py"


rule prepare_annot_pfamdb:
    input:  fa = "results/annotation/transdecoder_dir/longest_orfs.pep",
            pfam_db = os.path.join(GLOBAL_TMPD_PATH,"Pfam-A.hmm"),
    output: pfam_out = "results/annotation/TrinotatePFAM.out",
    log:    run = "logs/prepare_annot_pfamdb.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "results/annotation",
            hmm_out = "results/annotation/hmmscan.out",
    conda:  "../wrappers/prepare_annot_pfamdb/env.yaml"
    script: "../wrappers/prepare_annot_pfamdb/script.py"
    
    
rule prepare_annot_blastn:
    input:  fa = "results/assembly/merged/merged_samples.merged.fa",
            db = os.path.join(GLOBAL_REF_PATH,"general/nucleotide_DB/nt_db.ndb")
    output: blast_out = "results/annotation/blastn_{db}_processed.taxids_{tax}.tsv",
    log:    run = "logs/prepare_annot_blastn_to_{db}.taxids_{tax}.log"
    threads: 20
    resources:  mem = 50
    params: blast_full= "results/annotation/blastn_{db}_full.taxids_{tax}.tsv",
            eval_cutof = 1e-20, 
            pident_cutof = 80, # TODO: 0-100
            cov_cutof = 0.8,
            taxid = lambda wcs: [i.split(':')[1] for i in nt_taxids if i.split(':')[0].casefold() == wcs.tax.casefold()][0],
            tax_list = "results/annotation/taxids_{tax}",
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/prepare_annot_blastn/env.yaml"
    script: "../wrappers/prepare_annot_blastn/script.py"
    

def prepare_db_for_blast_rules(wcs):
    if "merops" in wcs.db:
        return os.path.join(GLOBAL_REF_PATH,"general/MEROPS_DB/merops_pepunit.diamond_db.dmnd")
    if "refseq" in wcs.db:
        return os.path.join(GLOBAL_REF_PATH,"general/RefSeq_DB/NCBI_RefSeq_complete_nonredundant_protein.diamond_db.dmnd")
    else:
        return os.path.join(GLOBAL_REF_PATH,"general/UniProt_SwissProt_DB/uniprot_sprot.dmnd")

rule prepare_annot_blastp:
    input:  fa = "results/annotation/transdecoder_dir/longest_orfs.pep",
            db = prepare_db_for_blast_rules,
    output: blast_out = "results/annotation/blastp_{db}_processed.tsv",
    log:    run = "logs/prepare_annot_blastp_to_{db}.log"
    threads: 20
    resources:  mem = 50
    params: blast_full= "results/annotation/blastp_{db}_full.tsv",
            eval_cutof = 1e-20,
            pident_cutof = 80, # 0-100
            cov_cutof = 0.8,
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/prepare_annot_blastp/env.yaml"
    script: "../wrappers/prepare_annot_blastp/script.py"


rule prepare_annot_blastx:
    input:  fa = "results/assembly/merged/merged_samples.merged.fa",
            db = prepare_db_for_blast_rules,
    output: blast_out = "results/annotation/blastx_{db}_processed.tsv",
    log:    run = "logs/prepare_annot_blastx_to_{db}.log"
    threads: 20
    resources:  mem = 50
    params: blast_full= "results/annotation/blastx_{db}_full.tsv",
            eval_cutof = 1e-20,
            pident_cutof = 80, # 0-100
            cov_cutof = 0.8,
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/prepare_annot_blastx/env.yaml"
    script: "../wrappers/prepare_annot_blastx/script.py"
    
    
# def prepare_annot_makeblastdb_input(wcs):
#     if "merops" in wcs.db:
#         return os.path.join(GLOBAL_REF_PATH,"/general/MEROPS_DB/merops_pepunit.fa.gz")
#     if "refseq" in wcs.db:
#         return os.path.join(GLOBAL_REF_PATH,"/general/RefSeq_DB/NCBI_RefSeq_complete_nonredundant_protein.fa.gz")
#     else:
#         return os.path.join(GLOBAL_REF_PATH,"/general/UniProt_SwissProt_DB/uniprot_sprot.fa.gz")
# 
# rule prepare_annot_makeblastdb(wcs):
#     input:  fa = prepare_annot_makeblastdb_input,
#     output: db_ok = os.path.join(GLOBAL_TMPD_PATH,"blastdb_{db}.ok"),
#     log:    run = "logs/prepare_annot_makeblastdb_for_{db}.log"
#     threads: 10
#     resources:  mem = 30
#     params: simple_names = os.path.join(GLOBAL_TMPD_PATH,"uniprot_sprot.simple_names.fasta"),
#             working_dir = "./",
#     conda:  "../wrappers/prepare_annot_makeblastdb/env.yaml"
#     script: "../wrappers/prepare_annot_makeblastdb/script.py"
    
    
rule prepare_annot_long_orfs:
    input:  fa = "results/assembly/merged/merged_samples.merged.fa",
            gene2trans = "results/assembly/merged/merged_samples.merged.fa.gene_trans_map",
    output: orfs_pep = "results/annotation/transdecoder_dir/longest_orfs.pep",
    log:    run = "logs/prepare_annot_long_orfs.log"
    threads: 1
    resources:  mem = 30
    params: stranded = config["trinity_strandness"],
            prefix = "results/annotation/",
    conda:  "../wrappers/prepare_annot_long_orfs/env.yaml"
    script: "../wrappers/prepare_annot_long_orfs/script.py"


rule prepare_annot_build:
    output: pfam_db = os.path.join(GLOBAL_TMPD_PATH,"Pfam-A.hmm"),
            sqlite = os.path.join(GLOBAL_TMPD_PATH,"Trinotate.sqlite"),
    log:    run = "logs/prepare_annot_build.log"
    threads: 1
    resources:  mem = 5
    params: prefix = GLOBAL_TMPD_PATH+"/",
            pfam_db = os.path.join(GLOBAL_TMPD_PATH,"Pfam-A.hmm.gz"),
    conda:  "../wrappers/prepare_annot_build/env.yaml"
    script: "../wrappers/prepare_annot_build/script.py"
