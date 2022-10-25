rule dummy:
    input:  expand("preprocess/mapped/{sample}.bam", sample = sample_tab.sample_name)
    output: "dummy"
    shell:  ""


def mark_duplicates_input(wildcards):
    input = {}
    input["bam"] = "preprocess/mapped/{sample}.not_markDups.bam"
    input["bai"] = "preprocess/mapped/{sample}.not_markDups.bam.bai"
    # input["transcriptome_bam"] = "preprocess/mapped/transcriptome/{sample}.not_markDups.transcriptome.bam"
    return input

rule mark_duplicates:
    input:  unpack(mark_duplicates_input)
    output: bam = "preprocess/mapped/{sample}.bam",
            bai = "preprocess/mapped/{sample}.bam.bai",
            # transcriptome_bam = "preprocess/mapped/transcriptome/{sample}.transcriptome.bam",
    log:    run = "logs/{sample}/mark_duplicates.log"
    threads: 8
    resources:  mem = 15
    params: mtx = "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
            mark_duplicates = config["mark_duplicates"],
            rmDup = config["remove_duplicates"], # allow possibility for rm duplicates true
            UMI = config["UMI"],
            umi_usage = config["umi_usage"],
            keep_not_markDups_bam = config["keep_not_markDups_bam"],
    conda:  "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


def align_to_genome_input(wildcards):
    preprocessed = "preprocess/trimmed_fastq"
    if read_pair_tags == ["SE"]:
        return os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]

rule align_to_genome:
    input:  fastqs = align_to_genome_input,
            genome = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"]),
            fai_ucsc = expand("{ref_dir}/seq/{ref}.fa.fai.ucsc",ref_dir=reference_directory,ref=config["reference"]),
            gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"]),
            index = expand("{ref_dir}/index/STAR/SAindex",ref_dir=reference_directory,ref=config["reference"])
    output: bam = "preprocess/mapped/{sample}.not_markDups.bam",
            bai = "preprocess/mapped/{sample}.not_markDups.bam.bai",
            # transcriptome_bam = "preprocess/mapped/transcriptome/{sample}.not_markDups.transcriptome.bam",
    log:    run = "logs/{sample}/align_to_genome.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "preprocess/mapped/{sample}/{sample}",
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
            tmpd = GLOBAL_TMPD_PATH,
    conda:  "../wrappers/align_to_genome/env.yaml"
    script: "../wrappers/align_to_genome/script.py"


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
            trim_left1 = config["trim_left1"], # Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
            trim_right1 = config["trim_right1"], # Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            trim_left2 = config["trim_left2"], # Applied only if trim left is true, trimming from R2 (different for classic:0, quant:?, sense:7)
            trim_right2 = config["trim_right2"], # Applied only if trim right is true, trimming from R2; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            phred = "-phred33",
            leading = 0,
            trailing = 0,
            crop = 250,
            minlen = config["min_length"],
            slid_w_1 = 4,
            slid_w_2 = 5,
    conda:  "../wrappers/trim_reads/env.yaml"
    script: "../wrappers/trim_reads/script.py"
