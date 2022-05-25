# This snakefile records the orpheum rules that were run on each spacegraphcats neighborhood.
# Orpheum predicts the open reading from of short sequencing reads and translates the read into protein space.
# Originally, we had planned to do metapangenome analysis (see https://github.com/taylorreiter/2021-paper-metapangenomes/).
# However, the goal of this analysis was to determine if sequences that were more abundant in IBD (CD or UC) were exclusive to IBD.
# I wanted to do this in protein space because it is more permissive to small evolutionary difference (e.g. third base pair wobble). 
# This ended up being difficult though, because we ran orpheum on the spacegraphcats neighborhoods from each metagenome, and ran differential abundance analysis on the metapangenome species graph (graph that contained sequences from all metagenomes, for one particular species).
# Unlike the metagenome queries, it wasn't as clear which specific reads led to each dominating set piece in the metapangenome species graph, and therefore it was hard to get the reads back that underlied the differentially abundant dominating set pieces. 
# Without reads, I couldn't run orpheum.
# What I ended up doing instead was using nucleotide k-mers from the dominating set pieces, because I didn't need to go back to the reads to get the translation.
# This ended up serving my purpose fine.

# This snakefile will only run if the main Snakefile is run.
# It relies on some output files from that pipeline; this is not a fully contained and automated pipeline, but only documents code that was run at some point.

import pandas as pd
#import feather
import numpy as np
from sourmash import signature
import sourmash
import glob
import os
import csv
import re
from collections import Counter

TMPDIR = "/scratch/tereiter/"

SEED = [1, 2, 3, 4, 5, 6]
ABUNDANCE = ['increased', 'decreased']

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
SAMPLES = m.sort_values(by='read_count')['run_accession']
LIBRARIES = m['library_name'].unique().tolist()
STUDY = m['study_accession'].unique().tolist()

class Checkpoint_AccToDbs:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_acc_dbs(self):
        acc_db_csv = f'outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.species_dbs.csv'
        assert os.path.exists(acc_db_csv)

        acc_dbs = []
        with open(acc_db_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['accession']
               db = row['species']
               acc_db = acc + "--"  + db
               acc_dbs.append(acc_db)

        return acc_dbs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'query_to_species_db';
        # this will trigger exception until that rule has been run.
        checkpoints.acc_to_species_db.get(**w)

        # parse accessions in gather output file
        genome_acc_dbs = self.get_acc_dbs()

        p = expand(self.pattern, acc_db=genome_acc_dbs, **w)
        return p

rule all:
    input:
        #Checkpoint_AccToDbs("outputs/sgc_genome_queries_orpheum_species_sketch_table/{acc_db}_long.csv"),
        #Checkpoint_AccToDbs("outputs/sgc_genome_queries_orpheum_species_comp/{acc_db}_clustered.csv"),

####################################################
## Metapangenome analysis of sgc genome query nbhds
####################################################

checkpoint acc_to_species_db:
    input: lineages = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.lineages.csv"
    output: csv = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.species_dbs.csv"
    conda: "envs/tidy.yml"
    resources:
        tmpdir = TMPDIR,
        mem_mb = 32000
    threads: 1
    script: "scripts/generate_acc_to_species_db.R" 
  
rule orpheum_translate_sgc_genome_query_nbhds:
    input: 
        ref = "inputs/orpheum_index/gtdb-rs202.{db}.protein-k10.nodegraph",
        fasta="outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    output:
        pep="outputs/sgc_genome_queries_orpheum_species/{library}-{acc}--{db}.coding.faa",
        nuc="outputs/sgc_genome_queries_orpheum_species/{library}-{acc}--{db}.nuc_coding.fna",
        nuc_noncoding="outputs/sgc_genome_queries_orpheum_species/{library}-{acc}--{db}.nuc_noncoding.fna",
        json="outputs/sgc_genome_queries_orpheum_species/{library}-{acc}--{db}.summary.json"
    conda: "envs/orpheum.yml"
    resources:  
        mem_mb=lambda wildcards, attempt: attempt * 5000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    orpheum translate --jaccard-threshold 0.39 --alphabet protein --peptide-ksize 10  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --json-summary {output.json} {input.ref} {input.fasta} > {output.pep}
    '''

rule sourmash_sketch_sgc_genome_query_nbhds_translated:
    input:"outputs/sgc_genome_queries_orpheum_species/{library}-{acc_db}.coding.faa"
    output: 'outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{acc_db}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.acc_db} {input}
    """

rule convert_signature_to_csv_species:
    input: 'outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{acc_db}.sig'
    output: 'outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{acc_db}.csv'
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_table_long_species:
    input: 
        expand('outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{{acc_db}}.csv', library = LIBRARIES)
    output: csv = "outputs/sgc_genome_queries_orpheum_species_sketch_table/{acc_db}_long.csv"
    conda: 'envs/tidy.yml'
    threads: 1
    resources:
        mem_mb=64000,
        tmpdir = TMPDIR
    script: "scripts/sketch_csv_to_long.R"

rule rename_signatures_species:
    input: 'outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{acc_db}.sig'
    output:'outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{acc_db}_renamed.sig'
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    shell:'''
    sourmash sig rename -o {output} {input} {wildcards.library}
    '''

rule compare_signatures_species:
    input: expand('outputs/sgc_genome_queries_orpheum_species_sigs/{library}-{{acc_db}}_renamed.sig', library = LIBRARIES) 
    output:
        comp="outputs/sgc_genome_queries_orpheum_species_comp/{acc_db}_comp",
        csv="outputs/sgc_genome_queries_orpheum_species_comp/{acc_db}_comp.csv"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=8000,
        tmpdir = TMPDIR
    shell:'''
    sourmash compare -o {output.comp} --csv {output.csv} {input}
    '''

rule plot_compare_signatures_species:
    input: "outputs/sgc_genome_queries_orpheum_species_comp/{acc_db}_comp",
    output: csv= "outputs/sgc_genome_queries_orpheum_species_comp/{acc_db}_clustered.csv"
    conda: 'envs/sourmash.yml'
    threads: 1
    params: outdir = "outputs/sgc_genome_queries_orpheum_species_comp/"
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    shell:'''
    sourmash plot --labels --pdf {input} --output-dir {params.outdir} --csv {output.csv}
    '''
