ruleorder: prepare_gene2trans > trinity_assembly

def multiqc_report_inputs(wc):
    inputs = {}
    inputs['rnaquast'] = "results/assembly/qc/rnaQUAST/short_report.txt"
    inputs['transrate']= "results/assembly/qc/transrate/assemblies.csv"
    inputs['stats'] = expand("results/assembly/trinity/{file_name}/merged_samples.{file_name}.general_stats.txt",
                              file_name = ['trinity_'+i for i in trinity_kmers])
    inputs['stats']+= expand("results/assembly/{file_name}/merged_samples.{file_name}.general_stats.txt",
                              file_name = ['merged'])
    inputs['busco'] = expand("results/assembly/qc/busco/{file_name}.{lineage}/{file_name}.short_summary.{lineage}.txt",
                            lineage = busco_filters.keys(),
                            file_name = ['merged']+['trinity_'+i for i in trinity_kmers])
    inputs['bowtie'] = expand("results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.flagstat",
                               file_name = ['merged']+['trinity_'+i for i in trinity_kmers])
    if config['genome_guided'] == False:
        inputs['busco']+= expand("results/assembly/qc/busco/{file_name}.{lineage}/{file_name}.short_summary.{lineage}.txt",
                                lineage = busco_filters.keys(),
                                file_name = ['megahit']+['spades_'+i for i in spades_kmers]+['trinity_'+i for i in trinity_kmers])
        inputs['bowtie']+= expand("results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.flagstat",
                                  file_name = ['megahit']+['spades_'+i for i in spades_kmers]+['trinity_'+i for i in trinity_kmers])
        inputs['stats']+= expand("results/assembly/{file_name}/merged_samples.{file_name}.general_stats.txt",
                                  file_name = ['megahit'])
        inputs['stats']+= expand("results/assembly/spades/{file_name}/merged_samples.{file_name}.general_stats.txt",
                                  file_name = ['spades_'+i for i in spades_kmers])
    return inputs

rule multiqc_report:
    input:  unpack(multiqc_report_inputs)
    output: html = "results/multiqc_report.html",
    log:    run = "logs/merged_samples/multiqc_report.log"
    threads:1
    params: multiqc_configs = [workflow.basedir+"/wrappers/multiqc_report/multiqc_config.yaml",\
                               workflow.basedir+"/wrappers/qc_assembly_rnaquast/multiqc_config.yaml",\
                               workflow.basedir+"/wrappers/qc_assembly_transrate/multiqc_config.yaml"],
            path = ".",
    conda:  "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"


rule qc_assembly_stats:
    input:  fa = "results/assembly/{file_name}.fa",
            gene2trans = "results/assembly/{file_name}.fa.gene_trans_map",
    output: trin_stats = "results/assembly/{file_name}.general_stats.txt",
    log:    run = "logs/merged_samples/{file_name}.qc_assembly_stats.log"
    threads: 1
    conda:  "../wrappers/qc_assembly_stats/env.yaml"
    script: "../wrappers/qc_assembly_stats/script.py"


rule prepare_gene2trans:
    input:  fa = "results/assembly/{filename}.fa",
    output: gene2trans = "results/assembly/{filename}.fa.gene_trans_map",
    log:    run = "logs/merged_samples/{filename}.prepare_annot_gene2trans.log"
    threads: 1
    resources:  mem = 10
    conda:  "../wrappers/prepare_gene2trans/env.yaml"
    script: "../wrappers/prepare_gene2trans/script.py"


def qc_assembly_input(ws):
    if "merged" in ws.file_name:
        return f"results/assembly/merged/merged_samples.{ws.file_name}.fa"
    elif "spades_" in ws.file_name:
        return f"results/assembly/spades/{ws.file_name}/merged_samples.{ws.file_name}.fa"
    elif "trinity_" in ws.file_name:
        return f"results/assembly/trinity/{ws.file_name}/merged_samples.{ws.file_name}.fa"
    elif "megahit" in ws.file_name:
        return f"results/assembly/megahit/merged_samples.{ws.file_name}.fa"
    else:
        return "ERROR"


rule qc_assembly_transrate:
    input:  merged = "results/assembly/merged/merged_samples.merged.fa",
            trinity = expand("results/assembly/trinity/trinity_{kmer}/merged_samples.trinity_{kmer}.fa", kmer = trinity_kmers),
            spades  = expand("results/assembly/spades/spades_{kmer}/merged_samples.spades_{kmer}.fa", kmer = spades_kmers),
            megahit = expand("results/assembly/megahit/merged_samples.megahit.fa"),
            r1 = "preprocess/corrected_fastq/merged_samples_R1.fastq.gz",
            r2 = "preprocess/corrected_fastq/merged_samples_R2.fastq.gz",
    output: summ = "results/assembly/qc/transrate/assemblies.csv",
    log:    run = "logs/merged_samples/qc_assembly_transrate.log",
    threads: 20
    params: prefix = "results/assembly/qc/transrate/",
            wdir = os.getcwd(),
            sdir = os.getcwd(),
            img = "docker://arnaudmeng/transrate:1.0.3",
            ref = config["transrate_ref"],
    # conda:  "../wrappers/qc_assembly_transrate/env.yaml"
    script: "../wrappers/qc_assembly_transrate/script.py"


rule qc_assembly_busco:
    input:  fa   = qc_assembly_input,
            busco= lambda wc: busco_filters[wc.lineage],
    output: summ = "results/assembly/qc/busco/{file_name}.{lineage}/{file_name}.short_summary.{lineage}.txt",
            ftab = "results/assembly/qc/busco/{file_name}.{lineage}/{file_name}.full_table.{lineage}.tsv",
            miss = "results/assembly/qc/busco/{file_name}.{lineage}/{file_name}.missing_busco.{lineage}.tsv",
    log:    run = "logs/merged_samples/qc_assembly_busco.{lineage}.{file_name}.log",
    threads: 10
    resources:  mem = 10
    params: prefix = "results/assembly/qc/busco/{file_name}.{lineage}/",
            #busco_db = os.path.join(GLOBAL_REF_PATH,"general/busco/odb10/"),
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/qc_assembly_busco/env.yaml"
    script: "../wrappers/qc_assembly_busco/script.py"


def qc_assembly_bowtie_input(wc):
    inputs = dict()
    inputs['fa'] = qc_assembly_input(wc)
    preprocessed = "preprocess/corrected_fastq"
    if read_pair_tags == ["SE"]:
        inputs['r1'] = f"{preprocessed}/merged_samples.fastq.gz"
    else:
        inputs['r1'] = f"{preprocessed}/merged_samples_R1.fastq.gz"
        inputs['r2'] = f"{preprocessed}/merged_samples_R2.fastq.gz"
    return inputs

rule qc_assembly_bowtie:
    input:  unpack(qc_assembly_bowtie_input),
    output: flagstat = "results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.flagstat",
            idxstats = "results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.idxstats",
            stats    = "results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.stats",
            picard_txt = "results/assembly/qc/bowtie/{file_name}/{file_name}.insert_size.txt",
            picard_pdf = "results/assembly/qc/bowtie/{file_name}/{file_name}.insert_size.pdf",
            ss_analysis = "results/assembly/qc/bowtie/{file_name}/{file_name}.strand_spec.dat.vioplot.pdf",
    log:    run = "logs/merged_samples/qc_assembly_bowtie.{file_name}.log"
    threads:40
    resources:  mem = 50
    params: prefix      = "results/assembly/qc/bowtie/{file_name}/{file_name}",
            bow_bam     = "results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.bam",
            bow_bam_tmp = "results/assembly/qc/bowtie/{file_name}/{file_name}.bw2.sam",
    conda:  "../wrappers/qc_assembly_bowtie/env.yaml"
    script: "../wrappers/qc_assembly_bowtie/script.py"


def qc_assembly_rnaquast_input(wc):
    inputs = {}
    inputs['genemark'] = os.path.join(GLOBAL_REF_PATH,"general/GeneMarkS-T.tar")
    inputs['assemblies'] = expand("results/assembly/merged/merged_samples.merged.fa")
    inputs['assemblies']+= expand("results/assembly/trinity/trinity_{kmer}/merged_samples.trinity_{kmer}.fa", kmer = trinity_kmers)
    if config['genome_guided'] == True:
        pass
    else:
        inputs['assemblies'] += expand("results/assembly/spades/spades_{kmer}/merged_samples.spades_{kmer}.fa", kmer = spades_kmers)
        inputs['assemblies'] += expand("results/assembly/megahit/merged_samples.megahit.fa")
    return inputs

rule qc_assembly_rnaquast:
    input:  unpack(qc_assembly_rnaquast_input),
    output: report = "results/assembly/qc/rnaQUAST/short_report.txt",
            full_rep = "results/assembly/qc/rnaQUAST/full_report_t.tsv",
            nx_plot = "results/assembly/qc/rnaQUAST/comparison_output/Nx_all.png",
            len_plot= "results/assembly/qc/rnaQUAST/comparison_output/transcript_length_all.png",
    log:    run = "logs/merged_samples/qc_assembly_rnaquast.log",
    threads:1
    params: tmpd     = GLOBAL_TMPD_PATH+"/",
            odir     = "results/assembly/qc/rnaQUAST/",
            stranded = config["trinity_strandness"],
    conda:  "../wrappers/qc_assembly_rnaquast/env.yaml"
    script: "../wrappers/qc_assembly_rnaquast/script.py"


def merge_assemblies_input(wc):
    inputs = {}
    inputs['evigene'] = os.path.join(GLOBAL_REF_PATH,"general/evigene.tar")
    inputs['assemblies'] = expand("results/assembly/trinity/trinity_{kmer}/merged_samples.trinity_{kmer}.fa", kmer = trinity_kmers)
    if config['genome_guided'] == True:
        pass
    else:
        inputs['assemblies'] += expand("results/assembly/spades/spades_{kmer}/merged_samples.spades_{kmer}.fa", kmer = spades_kmers)
        inputs['assemblies'] += expand("results/assembly/megahit/intermediate_contigs/k{kmer}.contigs.fa", kmer = megahit_kmers.split(','))
    return inputs

rule merge_assemblies:
    input:  unpack(merge_assemblies_input),
    output: trclass = "results/assembly/merged/merged_samples.merged.fa"
    log:    run = "logs/merged_samples/merge_assemblies.log"
    threads:40
    resources:  mem = 50
    params: base = "results/assembly/merged/merged_samples",
            base_dir = "results/assembly/merged/",
            trclass = "results/assembly/merged/okayset/merged_samples.compact.okay.compact.fa",
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/merge_assemblies/env.yaml"
    script: "../wrappers/merge_assemblies/script.py"


rule megahit_assembly:
    input:  r1 = expand("preprocess/corrected_fastq/merged_samples{read_tags}.fastq.gz",read_tags=pair_tag),
    output: fasta = "results/assembly/megahit/merged_samples.megahit.fa",
            contigs = expand("results/assembly/megahit/intermediate_contigs/k{kmer}.contigs.fa", kmer = megahit_kmers.split(',')),
    log:    run = "logs/merged_samples/megahit_assembly.log"
    threads: 20
    resources:  mem = 50
    params: kmers = megahit_kmers,
            fasta = "results/assembly/megahit/merged_samples.contigs.fa",
            mem = "0.2",
            paired = paired,
            tmpd = GLOBAL_TMPD_PATH+"/",
            odir = "results/assembly/megahit/",
            prefix = "merged_samples",
    conda:  "../wrappers/megahit_assembly/env.yaml"
    script: "../wrappers/megahit_assembly/script.py"


rule spades_assembly:
    input:  r1 = expand("preprocess/corrected_fastq/merged_samples{read_tags}.fastq.gz",read_tags=pair_tag),
    output: fasta = "results/assembly/spades/spades_{kmer}/merged_samples.spades_{kmer}.fa",
    log:    run = "logs/merged_samples/spades_assembly.kmer_{kmer}.log",
    threads: 20
    resources:  mem = 100
    params: stranded = config["trinity_strandness"],
            odir = "results/assembly/spades/spades_{kmer}",
            paired = paired,
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/spades_assembly/env.yaml"
    script: "../wrappers/spades_assembly/script.py"


def trinity_assembly_input(wc):
    inputs = {}
    if config['genome_guided'] == True:
        inputs['bam'] = f"preprocess/mapped/merged_samples.bam"
    else:
        inputs['r1'] = expand("preprocess/corrected_fastq/merged_samples{read_tags}.fastq.gz",read_tags=pair_tag)
    return inputs

rule trinity_assembly:
    input:  unpack(trinity_assembly_input),
    output: fasta = "results/assembly/trinity/trinity_{kmer}/merged_samples.trinity_{kmer}.fa",
            gtm = "results/assembly/trinity/trinity_{kmer}/merged_samples.trinity_{kmer}.fa.gene_trans_map",
    log:    run = "logs/merged_samples/trinity_assembly.kmer_{kmer}.log"
    threads: 20
    resources:  mem = 50
    params: stranded = config["trinity_strandness"],
            min_contig_len = config["trinity_min_contig_len"],
            min_kmer_cov = config["trinity_min_kmer_cov"],
            max_intron = config['max_intron'],
            min_map_cov = config['trinity_min_map_cov'],
            jaccard_clip = config['trinity_jaccard_clip'],
            use_genome = config['genome_guided'],
            paired = paired,
            tmpd = GLOBAL_TMPD_PATH+"/",
    conda:  "../wrappers/trinity_assembly/env.yaml"
    script: "../wrappers/trinity_assembly/script.py"
