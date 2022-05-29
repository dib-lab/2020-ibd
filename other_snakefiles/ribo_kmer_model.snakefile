import pandas as pd
import feather
import numpy as np
from sourmash import signature
import sourmash
import glob
import os
from collections import Counter

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
SAMPLES = m.sort_values(by='read_count')['run_accession']
LIBRARIES = m['library_name'].unique().tolist()
STUDY = m['study_accession'].unique().tolist()

rule all:
    input:
        expand('outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_validation_acc.csv', study = STUDY)

#rule compute_signatures_singlem:
#    input: "outputs/sgc_genome_queries_singlem_reads/{library}_singlem_reads.fq"
#    output: "outputs/sgc_genome_queries_singlem_sigs/{library}_singlem_reads.sig"
#    conda: "envs/sourmash.yml"
#    shell:'''
#    sourmash compute -k 31 --scaled 2000 -o {output} --track-abundance {input} || touch {output} 
#    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/sgc_genome_queries_singlem_sigs/{library}_singlem_reads.sig"
    output: "outputs/sgc_genome_queries_singlem_kmer_csv/{library}_singlem_reads.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output} || touch {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/sgc_genome_queries_singlem_kmer_csv/{library}_singlem_reads.csv", library = LIBRARIES)
    output: csv = "outputs/sgc_genome_queries_singlem_kmer_hash_tables/normalized_abund_hashes_long.csv"
    conda: 'envs/r.yml'
    script: "scripts/normalized_hash_abund_long_singlem.R"

rule make_hash_abund_table_wide:
    input: "outputs/sgc_genome_queries_singlem_kmer_hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/sgc_genome_queries_singlem_kmer_hash_tables/normalized_abund_hashes_wide.feather"
    run:
        import pandas as pd
        import feather
        
        ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
        ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
        ibd_wide = ibd_wide.fillna(0)
        ibd_wide['sample'] = ibd_wide.index
        ibd_wide = ibd_wide.reset_index(drop=True)
        ibd_wide.columns = ibd_wide.columns.astype(str)
        ibd_wide.to_feather(str(output)) 

########################################
## Random forests & optimization
########################################

rule install_pomona:
    input: "outputs/sgc_genome_queries_singlem_kmer_hash_tables/normalized_abund_hashes_wide.feather"
    output:
        pomona = "outputs/sgc_genome_queries_singlem_kmer_vita_rf/pomona_install.txt"
    conda: 'envs/rf.yml'
    script: "scripts/install_pomona.R"

rule vita_var_sel_rf:
    input:
        info = "inputs/working_metadata.tsv", 
        feather = "outputs/sgc_genome_queries_singlem_kmer_hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/sgc_genome_queries_singlem_kmer_vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/sgc_genome_queries_singlem_kmer_vita_rf/{study}_vita_rf.RDS",
        vita_vars = "outputs/sgc_genome_queries_singlem_kmer_vita_rf/{study}_vita_vars.txt",
        ibd_filt = "outputs/sgc_genome_queries_singlem_kmer_vita_rf/{study}_ibd_filt.csv"
    params: 
        threads = 4,
        validation_study = "{study}"
    conda: 'envs/rf.yml'
    script: "scripts/singlem_kmer_vita_rf.R"

rule loo_validation:
    input: 
        ibd_filt = 'outputs/sgc_genome_queries_singlem_kmer_vita_rf/{study}_ibd_filt.csv',
        info = 'inputs/working_metadata.tsv',
        eval_model = 'scripts/function_evaluate_model.R',
        ggconfusion = 'scripts/ggplotConfusionMatrix.R'
    output: 
        recommended_pars = 'outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_rec_pars.tsv',
        optimal_rf = 'outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_optimal_rf.RDS',
        training_accuracy = 'outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_training_acc.csv',
        training_confusion = 'outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_training_confusion.pdf',
        validation_accuracy = 'outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_validation_acc.csv',
        validation_confusion = 'outputs/sgc_genome_queries_singlem_kmer_optimal_rf/{study}_validation_confusion.pdf'
    params:
        threads = 4,
        validation_study = "{study}"
    conda: 'envs/tuneranger.yml'
    script: "scripts/tune_rf.R"
