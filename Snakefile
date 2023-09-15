import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]
#GLOBAL_RES_PATH = "/mnt/shared/CFBioinformatics/resources_backup"
#GLOBAL_REF_PATH = "/mnt/shared/CFBioinformatics/references_backup"
# GLOBAL_RES_PATH = "/mnt/ssd/ssd_3/resources/"
# GLOBAL_REF_PATH = "/mnt/references/"
# GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

#### Reference processing ####
# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
if 'reference' in config:
    config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
    reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

##### Config processing #####

# if not 'transrate_ref' in config:
#     config['transrate_ref'] = 'termites_proteins_ncbi.fasta'

if not 'min_qual' in config:
   config['min_qual'] = 25
if not 'remove_rRNAs' in config:
    config['remove_rRNAs'] = True # Remove rRNA reads or not (checkbox: [True, False], default: True)
if not 'overrep_qcov' in config:
    config['overrep_qcov'] = '0.66' # BLAST query coverage of overrepresented reads against rRNA DB index (real: 0.0 - 1.0, default: 0.66)
if not 'overrep_mism' in config:
    config['overrep_mism'] = '2' # Max number of mismatches of overrepresented reads against rRNA DB index (integer, default: 2)
if not 'overrep_kmer' in config:
    config['overrep_kmer'] = '50' # Jellyfish counting kmers of this size (integer, default: 50)
if not 'overrep_low_lim' in config:
    config['overrep_low_lim'] = '0.01' # Min fraction of overrepresented reads from total reads to consider for removing (real: 0.0 - 1.0, default: 0.01)
if not 'read_corr_kmer_len' in config:
    config["read_corr_kmer_len"] = 23 # Rcorrector's k-mer length (Integer: <=32, default: 23)
if not 'read_corr_window_max' in config:
    config["read_corr_window_max"] = 4 # Rcorrector's maximum number of correction within k-bp window (Integer, default: 4)
if not 'read_corr_weak_kmer_perc' in config:
    config["read_corr_weak_kmer_perc"] = 0.95 # Rcorrector's proportion of kmers that are used to estimate weak kmer count threshold (Real: 0.0 - 1.0, default: 0.95)
if not 'read_corr_100bp_max' in config:
    config["read_corr_100bp_max"] = 8 # Rcorrector's maximum number of correction every 100bp (Integer, default: 8)
if not 'trim_individually' in config:
    config['trim_individually'] = False
if not 'trinity_kmers' in config:
    config['trinity_kmers'] = "21;25;31"
if not 'trinity_min_contig_len' in config:
    config['trinity_min_contig_len'] = 100
if not 'trinity_min_kmer_cov' in config:
    config['trinity_min_kmer_cov'] = 2
if not 'trinity_jaccard_clip' in config:
    config['trinity_jaccard_clip'] = False
if not 'trinity_min_map_cov' in config:
    config['trinity_min_map_cov'] = 1
if not 'spades_kmers' in config:
    config['spades_kmers'] = '29;39;49;57;69;77;89;99;107;119;127'
if not 'megahit_kmers' in config:
    config['megahit_kmers'] = '21;29;39;59;79;99;119;141'
if not 'busco_lineage_version' in config:
    config['busco_lineage_version'] = 'odb10' # BUSCO lineage version appended to the lineage names
if not 'busco_lineage' in config:
    config['busco_lineage'] = "arachnida;arthropoda;metazoa" # BUSCO lineage names to check assembly against similar organisms
if not 'optional_prot_db_to_use' in config:
    config['optional_prot_db_to_use'] = 'merops' #TODO: add as a pipeline parameter for optional protein DBs from a list (refseq, merops, etc), uniprot is not optional, 
if not 'nt_taxids' in config:
    config['nt_taxids'] = 'Palpimanidae:336672;Araneae:6893' #TODO: add as a pipeline parameter for optional blast run against NT DB for limited tax subtree (the order matters, closer/smaller taxonomic units should be specified before the superior tax. units)

if config["trinity_kmers"] == "":
    trinity_kmers = '25'
else:
    trinity_kmers = config["trinity_kmers"].split(';')
spades_kmers = config["spades_kmers"].split(';')
megahit_kmers = config['megahit_kmers'].replace(';',',')

## set up BUSCO filters
busco_ref = pd.read_csv(os.path.join(GLOBAL_REF_PATH,"reference_info/busco_odb10_reference.txt"), sep=":")
ref_list = config['busco_lineage'].replace(' ','').split(";")
busco_filters = dict()
for r in ref_list:
    if r not in busco_ref["key"].unique():
        print("## WARNING: " + r + " is not a proper BUSCO reference key. It must be a key value from "+os.path.join(GLOBAL_REF_PATH,"reference_info/busco_odb10_reference.txt"))
    else:
        key = r+'_'+config['busco_lineage_version']
        path = os.path.join(GLOBAL_REF_PATH, busco_ref.loc[busco_ref["key"] == r, "path"].values[0])
        print("## Found BUSCO file: "+path)
        busco_filters[key] = path

## set up BLAST DBs
custom_prot_dbs = config['optional_prot_db_to_use'].replace(' ','').split(';')
all_prot_dbs = ['uniprot']+custom_prot_dbs
print("## INFO: Using following protein DBs as evidence sources: "+", ".join(all_prot_dbs))
nt_taxids = config['nt_taxids'].replace(' ','').split(';')
print("## INFO: Using following taxonomy units to narrow down the NCBI's nucleotide database: "+", ".join(nt_taxids))

sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
if not config['trim_individually']:
    sample_tab['ind_trim_left1'] = config['trim_left1']
    sample_tab['ind_trim_left2'] = config['trim_left2']
    sample_tab['ind_trim_right1'] = config['trim_right1']
    sample_tab['ind_trim_right2'] = config['trim_right2']
print(sample_tab)

if not config["is_paired"]:
    read_pair_tags = ["SE"]
    pair_tag = [""]
    paired = "SE"
else:
    read_pair_tags = ["R1","R2"]
    pair_tag = ["_R1","_R2"]
    paired = "PE"

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "R1|R2|SE"

# ##### Target rules #####

rule all:
    input:  "results/annotation/merged_samples.trinotate_annotation.gtf",
            "results/multiqc_report.html"

##### Modules #####

include: "rules/preprocess_reads.smk"
include: "rules/assemble_reads.smk"
include: "rules/annotate_assembly.smk"
