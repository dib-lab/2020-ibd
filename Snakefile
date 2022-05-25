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

class Checkpoint_GatherResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self):
        gather_csv = f'outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv'
        assert os.path.exists(gather_csv)

        genome_accs = []
        with open(gather_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['name'].split(' ')[0]
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {gather_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'gather_gtdb_rep_to_shared_assemblies'; 
        # this will trigger exception until that rule has been run.
        checkpoints.gather_gtdb_rep_to_shared_assemblies.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs()

        p = expand(self.pattern, acc=genome_accs, **w)
        return p


rule all:
    input:
        # SOURMASH COMPARE OUTPUTS:
        "outputs/comp/study_plt_all_filt_jaccard.pdf",
        "outputs/comp/diagnosis_plt_all_filt_jaccard.pdf",
        "outputs/comp/study_plt_all_filt_cosine.pdf",
        "outputs/comp/diagnosis_plt_all_filt_cosine.pdf",
        # VARIABLE SELECTION & RF OUTPUTS:
        "outputs/filt_sig_hashes/count_total_hashes.txt",
        expand('outputs/optimal_rf_seed/{study}_optimal_rf_seed{seed}.RDS', study = STUDY, seed = SEED),
        # VARIABLE CHARACTERIZATION OUTPUTS:
        expand("outputs/gather/{study}_vita_vars_gtdb_seed{seed}.csv", study = STUDY, seed = SEED),
        # SPACEGRAPHCATS OUTPUTS:
        Checkpoint_GatherResults("outputs/sgc_pangenome_catlases_corncob/{acc}_sig_ccs.tsv"),
        "outputs/sgc_genome_queries_fastp/all_fastp.tsv",
        #Checkpoint_GatherResults(expand("outputs/sgc_pangenome_catlases_corncob_sequences/{{acc}}_CD_{abundance}_contigs_search_gtdb_genomic.tsv", abundance = ABUNDANCE)),
        #Checkpoint_GatherResults(expand("outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_comp/{{acc}}_CD_{abundance}_contig_comp.csv", abundance = ABUNDANCE)),
        Checkpoint_GatherResults(expand("outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs_scaled1000.sig", abundance = ABUNDANCE)), 
        Checkpoint_GatherResults(expand("outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_intersect_long/{{acc}}_CD_{abundance}.csv", abundance = ABUNDANCE)),
        # CHARACTERIZING RESULTS OUTPUTS
        #expand("outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_all.csv", study = STUDY, seed = SEED),
        #expand("outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_pangenome_nbhd_reads.csv", study = STUDY, seed = SEED),
        #Checkpoint_GatherResults("outputs/sgc_pangenome_gather/{acc}_gtdb.csv"),

########################################
## PREPROCESSING
########################################

rule download_fastq_files_R1:
    output: 
        r1="inputs/raw/{sample}_1.fastq.gz",
    threads: 1
    resources:
        mem_mb=1000
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastq_1 = row['fastq_ftp_1'].values
        fastq_1 = fastq_1[0]
        shell("wget -O {output.r1} {fastq_1}")


rule download_fastq_files_R2:
    output:
        r2="inputs/raw/{sample}_2.fastq.gz"
    threads: 1
    resources:
        mem_mb=1000
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastq_2 = row['fastq_ftp_2'].values
        fastq_2 = fastq_2[0]
        shell("wget -O {output.r2} {fastq_2}")


rule cat_libraries_R1:
    input: expand("inputs/raw/{sample}_1.fastq.gz", sample = SAMPLES)
    output: expand("inputs/cat/{library}_1.fastq.gz", library = LIBRARIES)
    threads: 1
    resources:
        mem_mb=4000
    run: 
        merge_df = m[['library_name','run_accession']]
        merge_df = copy.deepcopy(merge_df)
        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_1.fastq.gz")
        merge_dict = merge_df.groupby('library_name')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
        for library in merge_dict.keys():
            # merge SRR files
            to_merge = merge_dict[library]
            # Check if the merged file results from a single or multiple fastq files.
            # For n-to-1 merging, concatenate input files to produce the output file
            merge_nb = len(to_merge)
            if merge_nb > 1:
                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + library + "_1.fastq.gz"
            else:
                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + library + "_1.fastq.gz"
            os.system(cmd)
    
rule cat_libraries_R2:
    input: expand("inputs/raw/{sample}_2.fastq.gz", sample = SAMPLES)
    output: expand("inputs/cat/{library}_2.fastq.gz", library = LIBRARIES)
    threads: 1
    resources:
        mem_mb=4000
    run: 
        merge_df = m[['library_name','run_accession']]
        merge_df = copy.deepcopy(merge_df)
        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_2.fastq.gz")
        merge_dict = merge_df.groupby('library_name')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
        for library in merge_dict.keys():
            # merge SRR files
            to_merge = merge_dict[library]
            # Check if the merged file results from a single or multiple fastq files.
            # For n-to-1 merging, concatenate input files to produce the output file
            merge_nb = len(to_merge)
            if merge_nb > 1:
                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + library + "_2.fastq.gz"
            else:
                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + library + "_2.fastq.gz"
            os.system(cmd)
    
rule adapter_trim_files:
    input:
        r1 = "inputs/cat/{library}_1.fastq.gz",
        r2 = 'inputs/cat/{library}_2.fastq.gz',
        adapters = 'inputs/adapters2.fa'
    output:
        r1 = temp('outputs/trim/{library}_R1.trim.fq.gz'),
        r2 = temp('outputs/trim/{library}_R2.trim.fq.gz'),
        o1 = temp('outputs/trim/{library}_o1.trim.fq.gz'),
        o2 = temp('outputs/trim/{library}_o2.trim.fq.gz')
    conda: 'envs/env.yml'
    threads: 1
    resources:
        mem_mb=8000
    shell:'''
     trimmomatic PE {input.r1} {input.r2} \
             {output.r1} {output.o1} {output.r2} {output.o2} \
             ILLUMINACLIP:{input.adapters}:2:0:15 MINLEN:31  \
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2
    '''

rule cutadapt_files:
    input:
        r1 = 'outputs/trim/{library}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{library}_R2.trim.fq.gz',
    output:
        r1 = temp('outputs/cut/{library}_R1.cut.fq.gz'),
        r2 = temp('outputs/cut/{library}_R2.cut.fq.gz'),
    conda: 'envs/env2.yml'
    threads: 1
    resources:
        mem_mb=8000
    shell:'''
    cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o {output.r1} -p {output.r2} {input.r1} {input.r2}
    '''

rule fastqc:
    input:
        r1 = 'outputs/cut/{library}_R1.cut.fq.gz',
        r2 = 'outputs/cut/{library}_R2.cut.fq.gz',
    output:
        r1 = 'outputs/fastqc/{library}_R1.cut_fastqc.html',
        r2 = 'outputs/fastqc/{library}_R2.cut_fastqc.html'
    conda: 'envs/env2.yml'
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    fastqc -o outputs/fastqc {input} 
    '''
    
rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    output:
        r1 = 'outputs/bbduk/{library}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{library}_R2.nohost.fq.gz',
        human_r1='outputs/bbduk/{library}_R1.human.fq.gz',
        human_r2='outputs/bbduk/{library}_R2.human.fq.gz'
    input: 
        r1 = 'outputs/cut/{library}_R1.cut.fq.gz',
        r2 = 'outputs/cut/{library}_R2.cut.fq.gz',
        human='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    threads: 1
    resources:
        mem_mb=64000
    conda: 'envs/env.yml'
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{library}_R1.nohost.fq.gz',
        'outputs/bbduk/{library}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    conda: 'envs/env.yml'
    threads: 1
    resources:
        mem_mb=64000
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule fastp_trimmed_reads:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: "outputs/fastp_abundtrim/{library}.abundtrim.fastp.json"
    conda: "envs/fastp.yml"
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    fastp -i {input} --interleaved_in -j {output}
    '''

rule multiqc_fastp:
    input: expand("outputs/fastp_abundtrim/{library}.abundtrim.fastp.json", library = LIBRARIES)
    output: "outputs/fastp_abundtrim/multiqc_data/mqc_fastp_filtered_reads_plot_1.txt"
    params: 
        indir = "outputs/fastp_abundtrim",
        outdir = "outputs/fastp_abundtrim"
    conda: "envs/multiqc.yml"
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
    '''

rule count_kmers_ntcard:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: 
        fstat = "outputs/ntcard/{library}.fstat",
        freq = "outputs/ntcard/{library}.freq"
    conda: 'envs/ntcard.yml'
    threads: 4
    resources:
        mem_mb=4000
    shell:'''
    ntcard -k31 -c2000 -t 4 -o {output.freq} {input} &> {output.fstat}
    '''

rule format_ntcard_kmer_count:
    input: fstat = expand("outputs/ntcard/{library}.fstat", library = LIBRARIES)
    output: tsv = 'outputs/ntcard/all_kmer_count.tsv'
    conda: "envs/tidy.yml"
    threads: 1
    resources:
        mem_mb=4000
    script: "scripts/format_ntcard_kmer_count.R"

rule compute_signatures:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: "outputs/sigs/{library}.sig"
    conda: 'envs/env.yml'
    threads: 1
    resources:
        mem_mb=2000
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} {input}
    '''

########################################
## Filtering and formatting signatures
########################################

rule get_greater_than_1_filt_sigs:
    input: expand("outputs/sigs/{library}.sig", library = LIBRARIES) 
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    threads: 1
    resources:
        mem_mb=64000
    run:
        # Determine the number of hashes, the number of unique hashes, and the number of
        # hashes that occur once across 605 gut metagenomes. Calculated for a scaled of 2k. 
        # 9 million hashes is the current approximate upper limit with which to build a 
        # sample vs hash abundance table using my current methods.

        files = input

        all_mins = []
        for file in files:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[1]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes 
                all_mins += mins

        counts = Counter(all_mins) # tally the number of hashes

        # remove hashes that occur only once
        for hashes, cnts in counts.copy().items():
            if cnts < 2:
                counts.pop(hashes)

        # write out hashes to a text file
        with open(str(output), "w") as f:
            for key in counts:
                print(key, file=f)


rule calc_total_hashes_sigs:
    """
    Output "statistics" about signatures
    """
    input: expand("outputs/sigs/{library}.sig", library = LIBRARIES)
    output: "outputs/filt_sig_hashes/count_total_hashes.txt"
    threads: 1
    resources:
        mem_mb=32000
    run:
        files = input
        all_mins = set()
        for file in files:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[1]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes
                all_mins.update(mins)

        with open(str(output), "w") as f:
            print(len(all_mins), file=f)


rule convert_greater_than_1_hashes_to_sig:
    input: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name greater_than_one_count_hashes --filename {input} {input}
    '''

rule filter_signatures_to_greater_than_1_hashes:
    input:
        filt_sig = "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig",
        sigs = "outputs/sigs/{library}.sig"
    output: "outputs/filt_sigs/{library}_filt.sig"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_filtered_sigs:
    input: "outputs/filt_sigs/{library}_filt.sig"
    output: "outputs/filt_sigs_named/{library}_filt_named.sig"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.library}_filt
    '''

rule describe_filtered_sigs:
    input: expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES)
    output: "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    sourmash signature describe --csv {output} {input}
    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/filt_sigs_named/{library}_filt_named.sig"
    output: "outputs/filt_sigs_named_csv/{library}_filt_named.csv"
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=2000
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/filt_sigs_named_csv/{library}_filt_named.csv", library = LIBRARIES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'envs/r.yml'
    threads: 1
    resources:
        mem_mb=64000
    script: "scripts/normalized_hash_abund_long.R"

rule make_hash_abund_table_wide:
    input: "outputs/hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    threads: 1
    resources:
        mem_mb=300000
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
    input: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    output:
        pomona = "outputs/vita_rf_seed/pomona_install.txt"
    conda: 'envs/rf.yml'
    threads: 1
    resources:
        mem_mb=1000
    script: "scripts/install_pomona.R"

rule vita_var_sel_rf_seed:
    input:
        info = "inputs/working_metadata.tsv", 
        feather = "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/vita_rf_seed/pomona_install.txt"
    output:
        vita_rf = "outputs/vita_rf_seed/{study}_vita_rf_seed{seed}.RDS",
        vita_vars = "outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.txt",
        ibd_filt = "outputs/vita_rf_seed/{study}_ibd_filt_seed{seed}.csv"
    resources:
        mem_mb=256000,
        runtime=12960
    threads: 32
    params: 
        threads = 32,
        validation_study = "{study}"
    conda: 'envs/rf.yml'
    script: "scripts/vita_rf_seed.R"

rule loo_validation_seed:
    input: 
        ibd_filt = 'outputs/vita_rf_seed/{study}_ibd_filt_seed{seed}.csv',
        info = 'inputs/working_metadata.tsv',
        eval_model = 'scripts/function_evaluate_model.R',
        ggconfusion = 'scripts/ggplotConfusionMatrix.R'
    output: 
        recommended_pars = 'outputs/optimal_rf_seed/{study}_rec_pars_seed{seed}.tsv',
        optimal_rf = 'outputs/optimal_rf_seed/{study}_optimal_rf_seed{seed}.RDS',
        training_accuracy = 'outputs/optimal_rf_seed/{study}_training_acc_seed{seed}.csv',
        training_confusion = 'outputs/optimal_rf_seed/{study}_training_confusion_seed{seed}.pdf',
        validation_accuracy = 'outputs/optimal_rf_seed/{study}_validation_acc_seed{seed}.csv',
        validation_confusion = 'outputs/optimal_rf_seed/{study}_validation_confusion_seed{seed}.pdf'
    resources:
        mem_mb = 16000,
        runtime=2880
    threads: 20
    params:
        threads = 20,
        validation_study = "{study}"
    conda: 'envs/tuneranger.yml'
    script: "scripts/tune_rf_seed.R"


############################################
## Predictive hash characterization - gather
############################################

rule convert_vita_vars_to_sig:
    input: "outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.txt"
    output: "outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.sig"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000
    threads: 1
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name vita_vars --filename {input} {input}
    '''

rule gather_vita_vars_gtdb:
    input:
        sig="outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.sig",
        db1="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic-reps.k31.sbt.zip",
        db2="/home/irber/sourmash_databases/outputs/sbt/genbank-viral-x1e6-k31.sbt.zip",
        db3="/home/irber/sourmash_databases/outputs/sbt/genbank-fungi-x1e6-k31.sbt.zip",
        db4="/home/irber/sourmash_databases/outputs/sbt/genbank-protozoa-x1e6-k31.sbt.zip",
    output: 
        csv="outputs/gather/{study}_vita_vars_gtdb_seed{seed}.csv",
        matches="outputs/gather/{study}_vita_vars_gtdb_seed{seed}.matches",
        un="outputs/gather/{study}_vita_vars_gtdb_seed{seed}.un"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 128000
    threads: 1
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db1} {input.db2} {input.db3} {input.db4}
    '''

# make this a checkpoint to interact with class Checkpoint_GatherResults
checkpoint gather_gtdb_rep_to_shared_assemblies:
    input:  
        gather=expand("outputs/gather/{study}_vita_vars_gtdb_seed{seed}.csv", study = STUDY, seed = SEED),
        gather_matches=expand("outputs/gather/{study}_vita_vars_gtdb_seed{seed}.matches", study = STUDY, seed = SEED),
        varimp=expand("outputs/optimal_rf_seed/{study}_optimal_rf_seed{seed}.RDS", study = STUDY, seed = SEED)
    output: 
        gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        gather_all_shared = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies_all.gather.csv"
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/gather_gtdb_rep_to_shared_assemblies.R"

# use to make acc:species db for orpheum open reading frame prediction (see orpheum*snakefile)
rule generate_shared_assembly_lineages:
    input:
        gather = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        db_lineages = "/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv"
    output: lineages = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.lineages.csv"
    conda: "envs/tidy.yml"
    resources: 
        mem_mb = 8000
    threads: 1
    script: "scripts/generate_shared_assembly_lineages.R"

rule generate_shared_assembly_lineages_all:
    input:
        gather = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies_all.gather.csv", 
        db_lineages = "/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv"
    output: lineages = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies_all.gather.csv"
    conda: "envs/tidy.yml"
    resources: 
        mem_mb = 8000
    threads: 1
    script: "scripts/generate_shared_assembly_lineages.R"

#############################################
# Spacegraphcats Genome Queries
#############################################

# specifying make_sgc_conf as the target downloads the genomes of interest,
# circumventing a bug in genome grist that prevents using the target
# download gather genomes. The sgc conf file is a dummy file -- it will be
# written to outputs/sgc, but the conf file has the wrong catlas bases.
rule download_shared_assemblies:
    input: 
        gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        conf = "inputs/genome-grist-conf.yml"
    output: "outputs/genbank_genomes_shared_assemblies/{acc}_genomic.fna.gz"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
    mv genbank_genomes/ outputs/genbank_genomes_shared_assemblies
    '''

rule generate_charcoal_genome_list:
    input:  ancient(Checkpoint_GatherResults("outputs/genbank_genomes_shared_assemblies/{acc}_genomic.fna.gz"))
    output: "outputs/charcoal_conf/charcoal.genome-list.txt"
    threads: 1
    resources:
        mem_mb=500
    shell:'''
    ls outputs/genbank_genomes_shared_assemblies/*gz | xargs -n 1 basename > {output} 
    '''

rule charcoal_decontaminate_shared_assemblies:
    input:
        genomes = ancient(Checkpoint_GatherResults("outputs/genbank_genomes_shared_assemblies/{acc}_genomic.fna.gz")),
        genome_list = "outputs/charcoal_conf/charcoal.genome-list.txt",
        conf = "inputs/charcoal-conf.yml",
        #genomes = "genbank_genomes/{acc}_genomic.fna.gz",
        #genome_list = "outputs/charcoal_conf/{acc}.genome-list.txt",
        #conf = "outputs/charcoal_conf/{acc}-conf.yml",
        genome_lineages = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.lineages.csv",
        db="/group/ctbrowngrp/gtdb/databases/gtdb-rs202.genomic.k31.zip",
        db_lineages="/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv"
    #output: "outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz"
    output: 
        hitlist="outputs/charcoal/stage1_hitlist.csv",
        clean_finished="outputs/charcoal/clean_finished.txt"
    resources:
        mem_mb = 128000
    threads: 8
    conda: "envs/charcoal.yml"
    shell:'''
    python -m charcoal run {input.conf} -j {threads} clean --nolock --latency-wait 15 --rerun-incomplete
    touch {output.clean_finished}
    '''

rule touch_decontaminated_shared_assemblies:
    input: 
        "outputs/charcoal/stage1_hitlist.csv",
        "outputs/charcoal/clean_finished.txt"
    output: "outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz"
    resources:
        mem_mb = 500
    threads: 1
    shell:'''
    touch {output}
    '''
    
rule make_sgc_genome_query_conf_files:
    input:
        csv = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        queries = ancient(Checkpoint_GatherResults("outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz")),
    output:
        conf = "outputs/sgc_conf/{library}_k31_r1_conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        query_list = "\n- ".join(input.queries)
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.library}
input_sequences:
- outputs/abundtrim/{wildcards.library}.abundtrim.fq.gz
ksize: 31
radius: 1
paired_reads: true
search:
- {query_list}
""", file=fp)

rule spacegraphcats_shared_assemblies:
    input: 
        queries = ancient(Checkpoint_GatherResults("outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz")), 
        conf = ancient("outputs/sgc_conf/{library}_k31_r1_conf.yml"),
        reads = "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output:
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/results.csv"
        #"outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz",
        #"outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.contigs.sig"
    params: outdir = "outputs/sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 500000
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    '''


# while this worked in the version of snakemake that I used, recent versions will erase the output file if it exists prior to running the rule.
# as such, this rule and the previous one should be replaced with a checkpoint approach.
# see here for example: https://github.com/taylorreiter/2022-infant-mge/blob/main/Snakefile#L414
rule touch_spacegraphcats_shared_assemblies:
    input: 
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/results.csv",
    output: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    resources:
        mem_mb = 500
    threads: 1
    shell:'''
    ls {output}
    '''

rule fastp_spacegraphcats_shared_assemblies:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    output: "outputs/sgc_genome_queries_fastp/{library}/{acc}_fastp.json"
    conda: "envs/fastp.yml"
    threads: 1
    resources:
        mem_mb=4000,
        tmpdir=TMPDIR
    shell:'''
    fastp -i {input} --interleaved_in -j {output}
    '''

rule multiqc_fastp_spacegraphcats_shared_assemblies:
    input: Checkpoint_GatherResults("outputs/sgc_genome_queries_fastp/{{library}}/{acc}_fastp.json")
    output: "outputs/sgc_genome_queries_fastp/{library}/multiqc_data/multiqc_data.json"
    params: 
        indir = lambda wildcards: "outputs/sgc_genome_queries_fastp/" + wildcards.library,
        outdir = lambda wildcards: "outputs/sgc_genome_queries_fastp/" + wildcards.library
    conda: "envs/multiqc.yml"
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    multiqc {params.indir} -o {params.outdir} --force 
    '''

rule summarize_multiqc_fastp_spacegraphcats_shared_assemblies:
    input: expand("outputs/sgc_genome_queries_fastp/{library}/multiqc_data/multiqc_data.json", library = LIBRARIES),
    output: tsv="outputs/sgc_genome_queries_fastp/all_fastp.tsv"
    conda: "envs/tidymultiqc.yml"
    threads: 1
    resources:
        mem_mb=16000
    script: "scripts/fastp_tidymultiqc.R"

####################################################
## Build spacegraphcats pangenome CAtlases
####################################################

rule diginorm_spacegraphcats_shared_assemblies:
    input: expand("outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{{acc}}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz", library = LIBRARIES)
    output: "outputs/sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    resources:
        mem_mb = 164000
    threads: 1
    conda: "envs/env.yml"
    shell:'''
    zcat {input} | normalize-by-median.py -k 20 -C 20 -M 164e9 --gzip -o {output} -
    '''

rule hardtrim_spacegraphcats_shared_assemblies:
    input: "outputs/sgc_genome_queries_diginorm/{acc}.diginorm.fa.gz"
    output: "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz"
    resources:
        mem_mb = 24000
    threads: 1
    conda: "envs/env.yml"
    shell:'''
    trim-low-abund.py -C 4 -M 20e9 -k 31 {input} --gzip -o {output}
    '''

rule make_sgc_pangenome_conf_files:
    input:
        reads = "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
    output:
        conf = "outputs/sgc_conf/{acc}_r10_conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.acc}
input_sequences:
- {input.reads}
radius: 10
paired_reads: true
""", file=fp)

rule spacegraphcats_pangenome_catlas_build:
    input:
        reads = "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
        conf = "outputs/sgc_conf/{acc}_r10_conf.yml"
    output: 
        "outputs/sgc_pangenome_catlases/{acc}_k31/cdbg.gxt",
        "outputs/sgc_pangenome_catlases/{acc}_k31/bcalm.unitigs.db",
        "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
    resources: mem_mb = 300000
    conda: "envs/spacegraphcats.yml"
    params: outdir = "outputs/sgc_pangenome_catlases"
    shell:'''
    python -m spacegraphcats build {input.conf} --outdir={params.outdir} --rerun-incomplete --nolock
    '''

rule spacegraphcats_pangenome_catlas_build_with_checkpoints:
# this rule is handy to viz the catlas, but otherwise is probably unnecessary
    input:
        reads = "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
        cdbg = "outputs/sgc_pangenome_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
    output: "outputs/sgc_pangenome_catlases/{acc}_k31_r10/10_1.checkpoint"
    resources: mem_mb = 100000
    conda: "envs/spacegraphcats.yml"
    params: 
        outdir = "outputs/sgc_pangenome_catlases",
        radius = 10,
        cdbg_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10", 
    shell:'''
    python -Werror -m spacegraphcats.catlas.catlas {params.cdbg_dir} {params.catlas_dir} {params.radius}
    '''

#rule spacegraphcats_pangenome_catlas_dom_graph:
#    shell:'''
#    # need to generalize/implement 02_visualize_sgc.ipynb
#    '''


#rule make_sgc_pangenome_multifasta_conf_files:
#    input:
#        reads = "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
#        ref_genes = "outputs/roary/{acc}/pan_genome_reference.fa",
#        ref_sig = "outputs/roary/{acc}/pan_genome_reference.sig"
#    output:
#        conf = "outputs/sgc_conf/{acc}_r10_multifasta_conf.yml"
#    resources:
#        mem_mb = 500
#    threads: 1
#    run:
#        query_list = "\n- ".join(input.queries)
#        with open(output.conf, 'wt') as fp:
#           print(f"""\
#catlas_base: {wildcards.acc}
#input_sequences:
#- {input.reads}
#radius: 10
#paired_reads: true
#multifasta_reference:
#- {input.ref_genes}
#multifasta_scaled: 2000
#multifasta_query_sig: {input.ref_sig}
#""", file=fp)
#
## TR TODO: UPDATE ENV? 
#rule spacegraphcats_pangenome_catlas_multifasta_annotate:
#    input:
#        reads = "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
#        ref_genes = "outputs/roary/{acc}/pan_genome_reference.fa",
#        ref_sig = "outputs/roary/{acc}/pan_genome_reference.sig",
#        conf = "outputs/sgc_conf/{acc}_r10_multifasta_conf.yml",
#        catlas = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
#    output: 
#        "outputs/sgc_pangenome_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_annot.csv",
#        "outputs/sgc_pangenome_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_by_record.csv"
#    params: 
#        outdir = "outputs/sgc_pangenome_catlases/",
#    conda: "envs/spacegraphcats_multifasta.yml"
#    resources:
#        mem_mb = 32000
#    threads: 1
#    shell:'''
#    python -m spacegraphcats {input.conf} multifasta_query --nolock --outdir {params.outdir} --rerun-incomplete
#    '''

rule spacegraphcats_pangenome_catlas_cdbg_to_pieces_map:
    input:
        cdbg = "outputs/sgc_pangenome_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
    output: "outputs/sgc_pangenome_catlases/{acc}_k31_r10/cdbg_to_pieces.csv"
    conda: "envs/spacegraphcats2.yml"
    resources: 
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10", 
    shell:'''
    scripts/cdbg_to_pieces.py {params.cdbg_dir} {params.catlas_dir}
    '''

rule tmp_cp_sgc_nbhds_w_lib_prefix:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    output: temp("outputs/sgc_genome_queries_tmp/{acc}/{library}.reads.gz")
    resources: 
        mem_mb = 500,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    cp {input} {output}
    '''
    
# TR TODO: update env to PR 303, or update sgc latest if merged. Since dom_abund is checked out, this might work like this...
rule spacegraphcats_pangenome_catlas_estimate_abundances:
    input:
        cdbg = "outputs/sgc_pangenome_catlases/{acc}_k31/cdbg.gxt",
        catlas = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv",
        reads = expand("outputs/sgc_genome_queries_tmp/{{acc}}/{library}.reads.gz", library = LIBRARIES)
        #reads = expand("outputs/abundtrim/{library}.abundtrim.fq.gz", library = LIBRARIES)
    output: expand("outputs/sgc_pangenome_catlases/{{acc}}_k31_r10_abund/{library}.reads.gz.dom_abund.csv", library = LIBRARIES)
    conda: "envs/spacegraphcats_dom.yml"
    resources: 
        mem_mb = 10000,
        tmpdir = TMPDIR
    threads: 1
    params:
        cdbg_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31" ,
        catlas_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10", 
        out_dir = lambda wildcards: "outputs/sgc_pangenome_catlases/" + wildcards.acc + "_k31_r10_abund", 
    shell:'''
    /home/tereiter/github/spacegraphcats/scripts/count-dominator-abundance.py {params.cdbg_dir} {params.catlas_dir} --outdir {params.out_dir} {input.reads}
    '''

rule format_spacegraphcats_pangenome_catlas_abundances:
    input: 
        dom_abund = expand("outputs/sgc_pangenome_catlases/{{acc}}_k31_r10_abund/{library}.reads.gz.dom_abund.csv", library = LIBRARIES)
    output: 
        dom_abund="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/all_dom_abund.tsv",
        dom_info="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/dom_info.tsv",
        dom_abund_pruned="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/all_dom_abund_pruned.tsv"
    conda: "envs/tidy.yml"
    resources: 
        mem_mb = 200000,
        tmpdir = TMPDIR
    threads: 1
    script: "scripts/format_pangenome_catlas_dom_abund.R"

rule install_corncob:
    output: corncob = "outputs/sgc_pangenome_catlases_corncob/corncob_install.txt"
    resources: 
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    conda: 'envs/corncob.yml'
    script: "scripts/install_corncob.R"

rule corncob_for_dominating_set_differential_abund:
    input: 
        corncob="outputs/sgc_pangenome_catlases_corncob/corncob_install.txt",
        dom_abund_pruned="outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/all_dom_abund_pruned.tsv",
        ntcard="outputs/ntcard/all_kmer_count.tsv",
        info = "inputs/working_metadata.tsv"
    output: 
        all_ccs = "outputs/sgc_pangenome_catlases_corncob/{acc}_all_ccs.tsv",
        sig_ccs = "outputs/sgc_pangenome_catlases_corncob/{acc}_sig_ccs.tsv"
    resources: 
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    conda: 'envs/corncob.yml'
    script: "scripts/corncob_dda.R"

rule grab_differentially_abundant_cdbg_ids:
    input: 
        sig_ccs = "outputs/sgc_pangenome_catlases_corncob/{acc}_sig_ccs.tsv",
        cdbg_to_pieces = "outputs/sgc_pangenome_catlases/{acc}_k31_r10/cdbg_to_pieces.csv",
    output:
        dom_ids_cd_increase = "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_increased_cdbg_ids.tsv.gz",
        dom_ids_cd_decrease = "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_decreased_cdbg_ids.tsv.gz",
        dom_ids_uc_increase = "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_UC_increased_cdbg_ids.tsv.gz",
        dom_ids_uc_decrease = "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_UC_decreased_cdbg_ids.tsv.gz",
    params: outdir = "outputs/sgc_pangenome_catlases_corncob_sequences/"
    resources: 
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    conda: 'envs/corncob.yml'
    script: "scripts/grab_differentially_abundant_dom_ids.R"

rule extract_contig_sequences_sig_cdbg_ids:
    # only run on cd increase for now; parameterize later if important to have for other sets
    input:
        contigs_db = "outputs/sgc_pangenome_catlases/{acc}_k31/bcalm.unitigs.db",
        cdbg_nbhds = "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_cdbg_ids.tsv.gz"
    output: "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.fa"
    conda: "envs/spacegraphcats2.yml"
    resources:
        mem_mb = 32000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    python -m spacegraphcats.search.extract_contigs --contigs-db {input.contigs_db} {input.cdbg_nbhds} -o {output}
    '''
 
rule sketch_contig_sequences_sig_cdbg_ids:
    input: "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.fa"
    output: "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.sig"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sourmash sketch dna -p k=31,scaled=2000 -o {output} {input}
    '''
 
rule sketch_contig_sequences_sig_cdbg_ids:
    input: "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.fa"
    output: "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs_scaled1000.sig"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sourmash sketch dna -p k=31,scaled=1000 -o {output} {input}
    '''
 
rule search_contig_sequences_sig_cdbg_ids:
    input: 
        sig="outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.sig",
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic.k31.sbt.zip"
    output: "outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs_search_gtdb_genomic.tsv"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 24000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sourmash search --threshold = 0.01 --max-containment --scaled 2000 -o {output} {input.sig} {input.db}
    '''

rule compare_diff_abund_contig_sequence_sigs_against_query_nbhds:
# this rule compares the differentially abundant sequences against sequences in the query nbhds
# to determine whether differentially abundant sequences are exclusive to IBD (CD, UC), or whether
# they occur across diagnosis conditions.
    input:
        query_nbhd_sigs=expand("outputs/sgc_genome_queries_nbhd_read_sigs/{library}/{{acc}}.cdbg_ids.reads.sig", library = LIBRARIES),
        diff_abund_sig="outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.sig",
    output: 
        comp="outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_comp/{acc}_CD_{abundance}_contig_comp",
        csv="outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_comp/{acc}_CD_{abundance}_contig_comp.csv"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 6000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sourmash compare -k 31 -o {output.comp} --csv {output.csv} {input}
    '''

rule intersect_diff_abund_contig_sequence_sigs_against_query_nbhds:
# this rule intersects the differentially abundant sequences against sequences in the query nbhds
# to determine which differentially abundant sequences are exclusive to IBD (CD, UC)
    input:
        query_nbhd_sig="outputs/sgc_genome_queries_nbhd_read_sigs/{library}/{acc}.cdbg_ids.reads.sig",
        diff_abund_sig="outputs/sgc_pangenome_catlases_corncob_sequences/{acc}_CD_{abundance}_contigs.sig",
    output: "outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_intersect/{library}/{acc}_CD_{abundance}.sig",
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sourmash sig intersect -k 31 -o {output} {input}
    '''

rule intersect_diff_abund_contig_sequence_sigs_against_query_nbhds_to_csv:
    input: "outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_intersect/{library}/{acc}_CD_{abundance}.sig",
    output: "outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_intersect/{library}/{acc}_CD_{abundance}.csv",
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule intersect_diff_abund_contig_sequence_sigs_against_query_nbhds_to_long:
    input: expand("outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_intersect/{library}/{{acc}}_CD_{{abundance}}.csv", library = LIBRARIES)
    output:csv="outputs/sgc_genome_queries_vs_pangenome_corncob_sequences_intersect_long/{acc}_CD_{abundance}.csv"
    conda: 'envs/tidy.yml'
    threads: 1
    resources:
        mem_mb=16000,
        tmpdir = TMPDIR
    script: "scripts/sketch_csv_to_long_name_lib.R"
 
##############################################
## Pangenome signature/variable importance
##############################################

# use default contig sigs output by sgc to start. 

rule calc_sig_sgc_genome_query_nbhd_reads:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
    output: "outputs/sgc_genome_queries_nbhd_read_sigs/{library}/{acc}.cdbg_ids.reads.sig"
    params: name = lambda wildcards: wildcards.library + "_" + wildcards.acc
    conda: "envs/sourmash.yml"
    resources:
        tmpdir = TMPDIR,
        mem_mb = 2000
    threads: 1
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} --merge {params.name} {input}
    '''

rule merge_sgc_sigs_to_pangenome:
    input: expand("outputs/sgc_genome_queries_nbhd_read_sigs/{library}/{{acc}}.cdbg_ids.reads.sig", library = LIBRARIES)
    output: "outputs/sgc_pangenome_nbhd_read_sigs/{acc}.sig"
    conda: "envs/sourmash.yml"
    resources:
        tmpdir = TMPDIR,
        mem_mb = 16000
    threads: 1
    shell:'''
    sourmash signature merge --name {wildcards.acc}_pangenome -o {output} -k 31 {input}
    '''

rule compare_sgc_pangenome_sigs:
    input: Checkpoint_GatherResults("outputs/sgc_pangenome_nbhd_read_sigs/{acc}.sig")
    output: 
        comp="outputs/sgc_pangenome_nbhd_read_compare/pangenome_compare.comp",
        csv="outputs/sgc_pangenome_nbhd_read_compare/pangenome_compare.csv"
    conda: "envs/sourmash.yml"
    resources:
        tmpdir = TMPDIR,
        mem_mb = 16000
    threads: 1
    shell:'''
    sourmash compare -o {output.comp} --csv {output.csv} --ignore-abundance {input}
    '''

rule index_sgc_pangenome_sigs:
    input: Checkpoint_GatherResults("outputs/sgc_pangenome_nbhd_read_sigs/{acc}.sig") 
    output: "outputs/sgc_pangenome_gather/sgc_pangenome_nbhd_reads_merged.sbt.json"
    conda: "envs/sourmash.yml"
    resources:
        tmpdir = TMPDIR,
        mem_mb = 32000
    threads: 1
    shell:'''
    sourmash index -k 31 {output} {input}
    '''

rule gather_vita_vars_study_against_sgc_pangenome_sigs:
    input:
        sig="outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.sig",
        db="outputs/sgc_pangenome_gather/sgc_pangenome_nbhd_reads_merged.sbt.json"
    output: 
        csv="outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_pangenome_nbhd_reads.csv",
        matches="outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_pangenome_nbhd_reads.matches",
        un="outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_pangenome_nbhd_reads.un",
    conda: 'envs/sourmash.yml'
    resources:
        tmpdir = TMPDIR,
        mem_mb = 16000
    threads: 1
    shell:'''
    sourmash gather --threshold-bp 0 -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

rule gather_against_sgc_pangenome_sigs_plus_all_dbs:
    input:
        sig="outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.sig",
        db0="outputs/sgc_pangenome_gather/sgc_pangenome_nbhd_reads_merged.sbt.json",
        db1="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic-reps.k31.sbt.zip"
    output: 
        csv="outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_all.csv",
        matches="outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_all.matches",
        un="outputs/sgc_pangenome_gather/{study}_vita_vars_seed{seed}_all.un"
    conda: 'envs/sourmash.yml'
    resources:
        tmpdir = TMPDIR,
        mem_mb = 32000
    threads: 1
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --output-unassigned {output.un} --save-matches {output.matches} --scaled 2000 -k 31 {input.sig} {input.db0} {input.db1}
    '''

rule gather_sgc_nbhds_against_gtdb:
    input:
        sig="outputs/sgc_pangenome_nbhd_read_sigs/{acc}.sig",
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic-reps.k31.sbt.zip"
    output: 
        csv="outputs/sgc_pangenome_gather/{acc}_gtdb.csv",
        matches="outputs/sgc_pangenome_gather/{acc}_gtdb.matches",
        un="outputs/sgc_pangenome_gather/{acc}_gtdb.un"
    conda: 'envs/sourmash.yml'
    resources:
        tmpdir = TMPDIR,
        mem_mb = 32000
    threads: 1
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --output-unassigned {output.un} --save-matches {output.matches} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

########################################
## PCoA
########################################

rule compare_signatures_cosine:
    input: 
        expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES),
    output: "outputs/comp/all_filt_comp_cosine.csv"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 32000
    threads: 8
    shell:'''
    sourmash compare -k 31 -p 8 --csv {output} {input}
    '''

rule compare_signatures_jaccard:
    input: 
        expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES),
    output: "outputs/comp/all_filt_comp_jaccard.csv"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 32000
    threads: 8
    shell:'''
    sourmash compare --ignore-abundance -k 31 -p 8 --csv {output} {input}
    '''

rule permanova_jaccard:
    input: 
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = "inputs/working_metadata.tsv",
        sig_info = "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    output: 
        perm = "outputs/comp/all_filt_permanova_jaccard.csv"
    conda: "envs/vegan.yml"
    resources:
        mem_mb = 16000
    threads: 1
    script: "scripts/run_permanova.R"

rule permanova_cosine:
    input: 
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = "inputs/working_metadata.tsv",
        sig_info = "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    output: 
        perm = "outputs/comp/all_filt_permanova_cosine.csv"
    conda: "envs/vegan.yml"
    resources:
        mem_mb = 16000
    threads: 1
    script: "scripts/run_permanova.R"

rule plot_comp_jaccard:
    input:
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        study = "outputs/comp/study_plt_all_filt_jaccard.pdf",
        diagnosis = "outputs/comp/diagnosis_plt_all_filt_jaccard.pdf"
    conda: "envs/ggplot.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/plot_comp.R"

rule plot_comp_cosine:
    input:
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        study = "outputs/comp/study_plt_all_filt_cosine.pdf",
        diagnosis = "outputs/comp/diagnosis_plt_all_filt_cosine.pdf"
    conda: "envs/ggplot.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/plot_comp.R"

########################################################
## Figures
########################################################

rule install_complexupset:
    output: complexupset = "Rmd_figures/complexupset_installed.txt"
    conda: "envs/rmd.yml"
    resources: 
        mem_mb = 1000
    threads: 1
    script: "scripts/install_complexupset.R"

#rule make_figures:
#    input:
#        complexupset="Rmd_figures/complexupset_installed.txt",
#        acc = expand("outputs/optimal_rf/{study}_{tv}_acc.csv", study = STUDY, tv = ['training', 'validation']), 
#        optimal_rf = expand("outputs/optimal_rf/{study}_optimal_rf.RDS", study = STUDY),
#        vita_rf = expand("outputs/vita_rf/{study}_vita_vars.txt", study = STUDY),
#        sgc_pangenome_gather_only = expand("outputs/sgc_pangenome_gather/{study}_vita_vars_pangenome.csv", study = STUDY),
#        sgc_pangenome_gather_all_db = expand("outputs/sgc_pangenome_gather/{study}_vita_vars_all.csv", study = STUDY),
#        gather_all_db = expand("outputs/gather/{study}_vita_vars_all.csv", study = STUDY),
#        gather_genbank = expand("outputs/gather/{study}_vita_vars_genbank.csv", study = STUDY),
#        gather_refseq = expand("outputs/gather/{study}_vita_vars_refseq.csv", study = STUDY),
#        hash_to_gather_map_loso = "outputs/gather_matches_loso_hash_map/hash_to_genome_map_at_least_5_studies.csv",
#        hash_to_gather_map_pangenome = "outputs/sgc_pangenome_gather/hash_to_genome_map_at_least_5_studies_pangenome.csv",
#        lca = "inputs/at_least_5_studies_vita_vars_gather_all_lca.csv",
#        gather_pangenome_vita = "outputs/sgc_pangenome_gather/at_least_5_studies_vita_vars_pangenome.csv"
#    output: "snakemake_figure_rmd.html"
#    conda: "envs/rmd.yml"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    script: "snakemake_figure_rmd.Rmd"  


#rule actually_make_figures:
#    input:
#        complexupset="Rmd_figures/complexupset_installed.txt",
#        acc = expand("outputs/optimal_rf/{study}_{tv}_acc.csv", study = STUDY, tv = ['training', 'validation']), 
#        optimal_rf = expand("outputs/optimal_rf/{study}_optimal_rf.RDS", study = STUDY),
#        vita_rf = expand("outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.txt", study = STUDY, seed = SEED),
#        sgc_pangenome_gather_only = expand("outputs/sgc_pangenome_gather/{study}_vita_vars_pangenome.csv", study = STUDY),
#        sgc_pangenome_gather_all_db = expand("outputs/sgc_pangenome_gather/{study}_vita_vars_all.csv", study = STUDY),
#        gather_all_db = expand("outputs/gather/{study}_vita_vars_all.csv", study = STUDY),
#        gather_genbank = expand("outputs/gather/{study}_vita_vars_genbank.csv", study = STUDY),
#        gather_refseq = expand("outputs/gather/{study}_vita_vars_refseq.csv", study = STUDY),
#    output: "figures_rmd.html"
#    conda: "envs/rmd.yml"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    Rscript -e "rmarkdown::render('figures_rmd.Rmd')"
#    '''
