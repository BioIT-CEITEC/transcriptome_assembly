
def mark_duplicates_input(wildcards):
    input = {}
    input["bam"] = "preprocess/mapped/merged_samples.not_markDups.bam"
    input["bai"] = "preprocess/mapped/merged_samples.not_markDups.bam.bai"
    return input

rule mark_duplicates:
    input:  unpack(mark_duplicates_input)
    output: bam = "preprocess/mapped/merged_samples.bam",
            bai = "preprocess/mapped/merged_samples.bam.bai",
    log:    run = "logs/merged_samples/mark_duplicates.log",
            mtx = "preprocess/mapped/merged_samples.markDups_metrics.txt",
    threads: 8
    resources:  mem = 15
    params: mark_duplicates = config["mark_duplicates"],
            rmDup = config["remove_duplicates"], # allow possibility for rm duplicates true
            UMI = config["UMI"],
            umi_usage = config["umi_usage"],
            keep_not_markDups_bam = config["keep_not_markDups_bam"],
            umi_stats_prefix = "preprocess/mapped/merged_samples.umi",
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


def align_to_genome_input(wc):
    inputs = {
        'genome': expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"]),
        'fai_ucsc': expand("{ref_dir}/seq/{ref}.fa.fai.ucsc",ref_dir=reference_directory,ref=config["reference"]),
        'gtf': expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"]),
        'index': expand("{ref_dir}/index/STAR/SAindex",ref_dir=reference_directory,ref=config["reference"])
    }
    preprocessed = "preprocess/trimmed_fastq"
    if read_pair_tags == ["SE"]:
        inputs['fastqs'] = [f"{preprocessed}/merged_samples.fastq.gz"]
    else:
        inputs['fastqs'] = [f"{preprocessed}/merged_samples_R1.fastq.gz", f"{preprocessed}/merged_samples_R2.fastq.gz"]
    return inputs

rule align_to_genome:
    input: unpack(align_to_genome_input)
    output: bam = "preprocess/mapped/merged_samples.not_markDups.bam",
            bai = "preprocess/mapped/merged_samples.not_markDups.bam.bai",
    log:    run = "logs/merged_samples/align_to_genome.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "preprocess/mapped/merged_samples/merged_samples",
            strandness = config["strandness"],
            num_mismatch= 999,  # Maximum number of mismatches; set this to high number (999) to disable and to use only perc_mismatch
            perc_mismatch= config["perc_mismatch"],
            max_intron= config["max_intron"],# Default used by ENCODE is 1000000; to turn this off set 1
            max_mate_dist=1000000,# Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this
            read_len=100,# Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
            organism=config["organism"],
            map_perc= config["map_perc"],
            map_score=config["map_score"],
            paired = paired,
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/align_to_genome/env.yaml"
    script: "../wrappers/align_to_genome/script.py"


rule merge_samples:
    input:  reads = expand("{{path}}/{sample}{{read_tags}}.fastq.gz", sample=sample_tab.sample_name),
            html  = expand("{{path}}/{sample}{{read_tags}}.fastqc.html", sample=sample_tab.sample_name),
    output: reads = "{path}/merged_samples{read_tags}.fastq.gz",
    log:    run = "{path}/merged_samples{read_tags}.merge_samples.log"
    conda:  "../wrappers/merge_samples/env.yaml"
    script: "../wrappers/merge_samples/script.py"


def read_correction_input(wcs):
  if config['remove_rRNAs']:
    if os.path.getsize(checkpoints.process_fastqc_data.get(sample=wcs.sample).output[0]) == 0:
        return {
            'reads' : expand("preprocess/removed_rRNAs_fastq/{sample}{read_tags}.fastq.gz",read_tags=pair_tag, sample=wcs.sample),
            'fastqc' : expand("preprocess/removed_rRNAs_fastq/{sample}{read_tags}.fastqc.html",read_tags=pair_tag, sample=wcs.sample),
        }
    else:
        return {
            'reads' : expand("preprocess/cleaned_fastq/{sample}{read_tags}.fastq.gz",read_tags=pair_tag, sample=wcs.sample),
            'fastqc' : expand("preprocess/cleaned_fastq/{sample}{read_tags}.fastqc.html",read_tags=pair_tag, sample=wcs.sample),
        }
  else:
    return {
      'reads' : expand("preprocess/trimmed_fastq/{sample}{read_tags}.fastq.gz",read_tags=pair_tag, sample=wcs.sample),
      'fastqc' : expand("preprocess/trimmed_fastq/{sample}{read_tags}.fastqc.html",read_tags=pair_tag, sample=wcs.sample),
    }

rule read_correction:
    input:  unpack(read_correction_input)
    output: reads = expand("preprocess/corrected_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    log:    run = "logs/{sample}/read_correction.log",
    threads:    20
    resources:  mem = 50
    params: mer_counts = "preprocess/corrected_fastq/{sample}.kmer_counts.jf",
            mer_counts_dump = "preprocess/corrected_fastq/{sample}.kmer_counts.fa",
            kmer = config["read_corr_kmer_len"], #23, # k-mer length (<=32, default: 23)
            wk = config["read_corr_weak_kmer_perc"], #0.95, # the proportion of kmers that are used to estimate weak kmer count threshold (Real, default: 0.95)
            max_corr = config["read_corr_100bp_max"], #8, # the maximum number of correction every 100bp (Integer, default: 8)
            max_corr_k = config["read_corr_window_max"], #4, # the maximum number of correction within k-bp window (Integer, default: 4)
            paired = config["is_paired"],
    conda:  "../wrappers/read_correction/env.yaml"
    script: "../wrappers/read_correction/script.py"


rule remove_overrepresented:
    input:  reads = expand("preprocess/removed_rRNAs_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
            tsv = "preprocess/removed_rRNAs_fastq/{sample}.overrep.toremove.tsv",
    output: reads = expand("preprocess/cleaned_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    log:    run = "logs/{sample}/remove_overrepresented.log",
    threads:20
    params: paired = config["is_paired"],
    conda:  "../wrappers/remove_overrepresented/env.yaml"
    script: "../wrappers/remove_overrepresented/script.py"


checkpoint process_fastqc_data:
    input:  db = os.path.join(GLOBAL_REF_PATH,"general/rRNA_DBs/rRNA_DB.blast_idx.done"),
            reads = expand("preprocess/removed_rRNAs_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
            fastqc = expand("preprocess/removed_rRNAs_fastq/{{sample}}{read_tags}.fastqc.html",read_tags=pair_tag),
            arch = expand("preprocess/removed_rRNAs_fastq/{{sample}}{read_tags}_fastqc.zip",read_tags=pair_tag)[0],
    output: tsv = "preprocess/removed_rRNAs_fastq/{sample}.overrep.toremove.tsv",
    log:    run = "logs/{sample}/process_fastqc_data.log",
    threads:20
    params: db = os.path.join(GLOBAL_REF_PATH,"general/rRNA_DBs/rRNA_DB.blast_idx"),
            paired = config["is_paired"],
            f1 = "preprocess/removed_rRNAs_fastq/{sample}.overrep.fa",
            b1 = "preprocess/removed_rRNAs_fastq/{sample}.overrep.blast.tsv",
            t1 = "preprocess/removed_rRNAs_fastq/{sample}.overrep.jellyfish.tsv",
            rscript = workflow.basedir + "/wrappers/process_fastqc_data/filter_blast_table.R",
            kmer = config["overrep_kmer"], # 50,
            low_lim = config["overrep_low_lim"], # 0.01,
            qcov = config["overrep_qcov"], # 0.66
            mism = config["overrep_mism"], # 2
    conda:  "../wrappers/process_fastqc_data/env.yaml"
    script: "../wrappers/process_fastqc_data/script.py"


rule remove_rRNAs:
    input:  reads = expand("preprocess/trimmed_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
            fastqc = expand("preprocess/trimmed_fastq/{{sample}}{read_tags}.fastqc.html",read_tags=pair_tag),
            db = os.path.join(GLOBAL_REF_PATH,"general/rRNA_DBs/rRNA_DB.bw2_idx.done"),
    output: reads = expand("preprocess/removed_rRNAs_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    log:    run = "logs/{sample}/remove_rRNAs.log",
            meta = "preprocess/removed_rRNAs_fastq/{sample}.bowtie2_metrics.txt",
    threads:    20
    resources:  mem = 30
    params: bw2_idx = os.path.join(GLOBAL_REF_PATH,"general/rRNA_DBs/rRNA_DB.bw2_idx"),
            paired = config["is_paired"],
            pu = "preprocess/removed_rRNAs_fastq/{sample}_R%.fastq.gz",
            uu = "preprocess/removed_rRNAs_fastq/{sample}.fastq.gz",
            stranded = config["trinity_strandness"], # RF,
            phred = "--phred33",
    conda:  "../wrappers/remove_rRNAs/env.yaml"
    script: "../wrappers/remove_rRNAs/script.py"


rule run_fastqc:
    input:  fastq= "{path}/{filename}.fastq.gz"
    output: html = "{path}/{filename}.fastqc.html",
            arch = "{path}/{filename}_fastqc.zip",
    log:    run  = "{path}/{filename}.run_fastqc.log"
    params: extra = "--noextract --format fastq --nogroup",
            tmp = GLOBAL_TMPD_PATH+"/",
            html = "{path}/{filename}_fastqc.html"
    conda:  "../wrappers/run_fastqc/env.yaml"
    script: "../wrappers/run_fastqc/script.py"


rule trim_reads:
    input:  raw = expand("raw_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    output: cleaned = expand("preprocess/trimmed_fastq/{{sample}}{read_tags}.fastq.gz",read_tags=pair_tag),
    log:    run = "logs/{sample}/trim_reads.log",
            trim = "preprocess/trimmed_fastq/{sample}.trim_stats.log",
    threads: 10
    resources:  mem = 10
    params: adaptors = config["trim_adapters"],
            r1u = "preprocess/trimmed_fastq/{sample}_R1.singletons.fastq.gz",
            r2u = "preprocess/trimmed_fastq/{sample}_R2.singletons.fastq.gz",
            trim_left1  = lambda wcs: sample_tab.loc[sample_tab.sample_name == wcs.sample, "ind_trim_left1"].min(), # Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
            trim_right1 = lambda wcs: sample_tab.loc[sample_tab.sample_name == wcs.sample, "ind_trim_right1"].min(), # Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            trim_left2  = lambda wcs: sample_tab.loc[sample_tab.sample_name == wcs.sample, "ind_trim_left2"].min(), # Applied only if trim left is true, trimming from R2 (different for classic:0, quant:?, sense:7)
            trim_right2 = lambda wcs: sample_tab.loc[sample_tab.sample_name == wcs.sample, "ind_trim_right2"].min(), # Applied only if trim right is true, trimming from R2; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            phred = "-phred33",
            leading = config['min_qual'],
            trailing = config['min_qual'],
            crop = 250,
            minlen = config["min_length"],
            slid_w_1 = 5,
            slid_w_2 = config['min_qual'],
    conda:  "../wrappers/trim_reads/env.yaml"
    script: "../wrappers/trim_reads/script.py"
