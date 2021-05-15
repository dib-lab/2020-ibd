import pandas as pd
#import feather
import numpy as np
from sourmash import signature
import sourmash
import glob
import os
from collections import Counter

SEED = [1, 2, 3, 4, 5, 6]

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
SAMPLES = m.sort_values(by='read_count')['run_accession']
LIBRARIES = m['library_name'].unique().tolist()
STUDY = m['study_accession'].unique().tolist()

# this variable is output after random forests. 
# Specifying this variable can be avoided by using checkpoints,
# but this makes the DAG take forever to solve. 
# this step requires user input anyway to download the matching
# random forest genomes, so specifying this variable manually is a 
# compromise. 
GATHER_GENOMES = ["ERS235530_10.fna", "ERS235531_43.fna", "ERS235603_16.fna", 
           "XieH_2016__YSZC12003_37172__bin.63.fa", "ZeeviD_2015__PNP_Main_232__bin.27.fa"]

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
        #"outputs/gather_matches_hash_map/hash_to_genome_map_gather_all.csv",
        # SPACEGRAPHCATS OUTPUTS:
        expand("outputs/nbhd_reads_sigs_csv/{library}/{gather_genome}.cdbg_ids.reads.csv", library = LIBRARIES, gather_genome = GATHER_GENOMES),
        #expand("outputs/sgc_genome_queries/{library}_k31_r1_multifasta/query-results.csv", library = LIBRARIES),
        # PANGENOME SIGS
        #expand("outputs/sgc_pangenome_gather/{study}_vita_vars_all.csv", study = STUDY),
        #expand("outputs/sgc_pangenome_gather/{study}_vita_vars_pangenome.csv", study = STUDY),
        #"figures_rmd.html",
        #"outputs/gather_matches_loso_multifasta/all-multifasta-query-results.emapper.annotations",
        # SINGLEM OUTPUTS:
        expand('outputs/singlem_abundtrim_optimal_rf/{study}_validation_acc.csv', study = STUDY),
        #expand('outputs/singlem_optimal_rf/{study}_validation_acc.csv', study = STUDY),
        #expand('outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_validation_acc.csv', study = STUDY),
        expand('outputs/singlem_abundtrim_kmer_optimal_rf/{study}_validation_acc.csv', study = STUDY)

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
        r1 = 'outputs/trim/{library}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{library}_R2.trim.fq.gz',
        o1 = 'outputs/trim/{library}_o1.trim.fq.gz',
        o2 = 'outputs/trim/{library}_o2.trim.fq.gz'
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
        r1 = 'outputs/cut/{library}_R1.cut.fq.gz',
        r2 = 'outputs/cut/{library}_R2.cut.fq.gz',
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
## Try singlem on abundtrim
########################################

rule split_paired_reads_abundtrim:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: 
        O = "outputs/abundtrim_split/{library}_orphan.abundtrim.fq.gz",
        R1 = "outputs/abundtrim_split/{library}_R1.abundtrim.fq.gz",
        R2 = "outputs/abundtrim_split/{library}_R2.abundtrim.fq.gz"
    conda: "envs/env.yml"
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    split-paired-reads.py -0 {output.O} -1 {output.R1} -2 {output.R2} --gzip {input}    
    '''

rule singlem_default_abundtrim:
    input: 
        R1 = "outputs/abundtrim_split/{library}_R1.abundtrim.fq.gz",
        R2 = "outputs/abundtrim_split/{library}_R2.abundtrim.fq.gz"
    output: "outputs/singlem_abundtrim/{library}_otu_default.csv"
    conda: "envs/singlem.yml"
    threads: 2
    resources:
        mem_mb=4000
    shell: '''
    singlem pipe --forward {input.R1} --reverse {input.R2} --otu_table {output} --output_extras --threads {threads} --filter_minimum_nucleotide 36 --min_orf_length 36 --filter_minimum_protein 12 #|| touch {output}
    touch {output} # creates output file for runs with no seq matches
    '''

rule singlem_16s_abundtrim_R1:
    input: 
        R1 = "outputs/abundtrim_split/{library}_R1.abundtrim.fq.gz",
        pkg =  "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg/CONTENTS.json"
    output: "outputs/singlem_abundtrim/{library}_otu_16s_R1.csv"
    conda: "envs/singlem.yml"
    threads: 2
    resources:
        mem_mb=4000
    params: 
        pkg_dir = "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg"
    shell: '''
    singlem pipe --sequences {input.R1} --singlem_packages {params.pkg_dir} --otu_table {output} --output_extras --threads {threads} --filter_minimum_nucleotide 36 --min_orf_length 36 --filter_minimum_protein 12 # || touch {output}
    touch {output}
    '''

rule singlem_16s_abundtrim_R2:
    input: 
        R2 = "outputs/abundtrim_split/{library}_R2.abundtrim.fq.gz",
        pkg =  "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg/CONTENTS.json"
    output: "outputs/singlem_abundtrim/{library}_otu_16s_R2.csv"
    conda: "envs/singlem.yml"
    threads: 2
    resources:
        mem_mb=4000
    params: 
        pkg_dir = "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg"
    shell: '''
    singlem pipe --sequences {input.R2} --singlem_packages {params.pkg_dir} --otu_table {output} --output_extras --threads {threads} --filter_minimum_nucleotide 36 --min_orf_length 36 --filter_minimum_protein 12 # || touch {output}
    touch {output}
    '''

rule combine_singlem_abundtrim:
    input:
        default = expand("outputs/singlem_abundtrim/{library}_otu_default.csv", library = LIBRARIES),
        s16_R1 = expand("outputs/singlem_abundtrim/{library}_otu_16s_R1.csv", library = LIBRARIES),
        s16_R2 = expand("outputs/singlem_abundtrim/{library}_otu_16s_R2.csv", library = LIBRARIES)
    output: res = "outputs/singlem_abundtrim/combined.tsv"
    conda: "envs/tidy.yml"
    threads: 1
    resources:
        mem_mb=16000
    script: "scripts/parse_singlem_abundtrim.R"

rule singlem_to_counts_abuntrim:
    input:  res = "outputs/singlem_abundtrim/combined.tsv"
    output: counts = "outputs/singlem_abundtrim/singlem_counts.tsv"
    conda: "envs/tidy.yml"
    threads: 1
    resources:
        mem_mb=16000
    script: "scripts/make_singlem_counts_abundtrim.R"

rule singlem_abundtrim_install_pomona:
    input: "outputs/singlem_abundtrim/combined.tsv"
    output:
        pomona = "outputs/singlem_abundtrim_vita_rf/pomona_install.txt"
    conda: 'envs/rf.yml'
    threads: 1
    resources:
        mem_mb=1000
    script: "scripts/install_pomona.R"

rule singlem_abundtrim_var_sel_rf:
    input:
        info = "inputs/working_metadata.tsv", 
        counts = "outputs/singlem_abundtrim/singlem_counts.tsv",
        pomona = "outputs/singlem_abundtrim_vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/singlem_abundtrim_vita_rf/{study}_vita_rf.RDS",
        vita_vars = "outputs/singlem_abundtrim_vita_rf/{study}_vita_vars.txt",
        ibd_filt = "outputs/singlem_abundtrim_vita_rf/{study}_ibd_filt.csv"
    threads: 6
    resources:
        mem_mb=32000
    params: 
        threads = 6,
        validation_study = "{study}"
    conda: 'envs/rf.yml'
    script: "scripts/singlem_vita_rf.R"

rule singlem_abundtrim_loo_validation:
    input: 
        ibd_filt = 'outputs/singlem_abundtrim_vita_rf/{study}_ibd_filt.csv',
        info = 'inputs/working_metadata.tsv',
        eval_model = 'scripts/function_evaluate_model.R',
        ggconfusion = 'scripts/ggplotConfusionMatrix.R'
    output: 
        recommended_pars = 'outputs/singlem_abundtrim_optimal_rf/{study}_rec_pars.tsv',
        optimal_rf = 'outputs/singlem_abundtrim_optimal_rf/{study}_optimal_rf.RDS',
        training_accuracy = 'outputs/singlem_abundtrim_optimal_rf/{study}_training_acc.csv',
        training_confusion = 'outputs/singlem_abundtrim_optimal_rf/{study}_training_confusion.pdf',
        validation_accuracy = 'outputs/singlem_abundtrim_optimal_rf/{study}_validation_acc.csv',
        validation_confusion = 'outputs/singlem_abundtrim_optimal_rf/{study}_validation_confusion.pdf'
    threads: 6
    resources:
        mem_mb=8000
    params:
        threads = 6,
        validation_study = "{study}"
    conda: 'envs/tuneranger.yml'
    script: "scripts/singlem_tune_rf.R"

######## kmer model of singlem on abundtrim

######### default all

rule extract_singlem_read_names_default_abundtrim:
    input: singlem = "outputs/singlem_abundtrim/{library}_otu_default.csv",
    output: reads = "outputs/singlem_abundtrim_reads/{library}_otu_default_read_names.txt" 
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/extract_singlem_read_names.R"

rule extract_singlem_reads_default:
    input: 
        names = "outputs/singlem_abundtrim_reads/{library}_otu_default_read_names.txt",
        fq = "outputs/abundtrim/{library}.abundtrim.fq.gz",
    output: "outputs/singlem_abundtrim_reads/{library}_otu_default.fq",
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    scripts/extract-aaseq-matches.py {input.names} {input.fq} > {output}
    '''

####### 16s R1
rule extract_singlem_read_names_16s_R1_abundtrim:
    input: singlem = "outputs/singlem_abundtrim/{library}_otu_16s_R1.csv",
    output: reads = "outputs/singlem_abundtrim_reads/{library}_otu_16s_R1_read_names.txt" 
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/extract_singlem_read_names.R"

rule extract_singlem_reads_16s_R1_abundtrim:
    input: 
        names = "outputs/singlem_abundtrim_reads/{library}_otu_16s_R1_read_names.txt",
        fq = "outputs/abundtrim/{library}.abundtrim.fq.gz",
    output: "outputs/singlem_abundtrim_reads/{library}_otu_16s_R1.fq",
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    scripts/extract-aaseq-matches.py {input.names} {input.fq} > {output}
    '''

######### 16 R2

rule extract_singlem_read_names_16s_R2_abundtrim:
    input: singlem = "outputs/singlem_abundtrim/{library}_otu_16s_R2.csv",
    output: reads = "outputs/singlem_abundtrim_reads/{library}_otu_16s_R2_read_names.txt" 
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/extract_singlem_read_names.R"

rule extract_singlem_reads_16s_R2_abundtrim:
    input: 
        names = "outputs/singlem_abundtrim_reads/{library}_otu_16s_R2_read_names.txt",
        fq = "outputs/abundtrim/{library}.abundtrim.fq.gz",
    output: "outputs/singlem_abundtrim_reads/{library}_otu_16s_R2.fq",
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    scripts/extract-aaseq-matches.py {input.names} {input.fq} > {output}
    '''

rule combine_singlem_reads_per_lib_abundtrim:
    input:
        "outputs/singlem_abundtrim_reads/{library}_otu_default.fq",
        "outputs/singlem_abundtrim_reads/{library}_otu_16s_R1.fq",
        "outputs/singlem_abundtrim_reads/{library}_otu_16s_R2.fq",
    output: "outputs/singlem_abundtrim_reads/{library}_singlem_reads.fq"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    cat {input} > {output}
    '''

rule compute_signatures_singlem:
    input: "outputs/singlem_abundtrim_reads/{library}_singlem_reads.fq"
    output: "outputs/singlem_abundtrim_sigs/{library}_singlem_reads.sig"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000
    threads: 1
    shell:'''
    sourmash compute -k 31 --scaled 2000 -o {output} --track-abundance {input} || touch {output} 
    '''

rule convert_greater_than_1_signatures_to_csv_abundtrim:
    input: "outputs/singlem_abundtrim_sigs/{library}_singlem_reads.sig"
    output: "outputs/singlem_abundtrim_kmer_csv/{library}_singlem_reads.csv"
    conda: 'envs/sourmash.yml'
    resources: 
        mem_mb=1000
    threads: 1
    shell:'''
    python scripts/sig_to_csv.py {input} {output} || touch {output}
    '''

rule make_hash_abund_table_long_normalized_abundtrim:
    input: 
        expand("outputs/singlem_abundtrim_kmer_csv/{library}_singlem_reads.csv", library = LIBRARIES)
    output: csv = "outputs/singlem_abundtrim_kmer_hash_tables/normalized_abund_hashes_long.csv"
    conda: 'envs/r.yml'
    resources: 
        mem_mb=16000
    threads: 1
    script: "scripts/normalized_hash_abund_long_singlem.R"

rule make_hash_abund_table_wide_abundtrim:
    input: "outputs/singlem_abundtrim_kmer_hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/singlem_abundtrim_kmer_hash_tables/normalized_abund_hashes_wide.feather"
    resources: 
        mem_mb=32000
    threads: 1
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

rule singlem_kmer_install_pomona_abundtrim:
    input: "outputs/singlem_abundtrim_kmer_hash_tables/normalized_abund_hashes_wide.feather"
    output:
        pomona = "outputs/singlem_abundtrim_kmer_vita_rf/pomona_install.txt"
    conda: 'envs/rf.yml'
    resources: 
        mem_mb=1000
    threads: 1
    script: "scripts/install_pomona.R"

rule singlem_kmer_vita_var_sel_rf_abundtrim:
    input:
        info = "inputs/working_metadata.tsv", 
        feather = "outputs/singlem_abundtrim_kmer_hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/singlem_abundtrim_kmer_vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/singlem_abundtrim_kmer_vita_rf/{study}_vita_rf.RDS",
        vita_vars = "outputs/singlem_abundtrim_kmer_vita_rf/{study}_vita_vars.txt",
        ibd_filt = "outputs/singlem_abundtrim_kmer_vita_rf/{study}_ibd_filt.csv"
    resources: 
        mem_mb=16000
    threads: 4
    params: 
        threads = 4,
        validation_study = "{study}"
    conda: 'envs/rf.yml'
    script: "scripts/singlem_kmer_vita_rf.R"

rule singlem_kmer_loo_validation:
    input: 
        ibd_filt = 'outputs/singlem_abundtrim_kmer_vita_rf/{study}_ibd_filt.csv',
        info = 'inputs/working_metadata.tsv',
        eval_model = 'scripts/function_evaluate_model.R',
        ggconfusion = 'scripts/ggplotConfusionMatrix.R'
    output: 
        recommended_pars = 'outputs/singlem_abundtrim_kmer_optimal_rf/{study}_rec_pars.tsv',
        optimal_rf = 'outputs/singlem_abundtrim_kmer_optimal_rf/{study}_optimal_rf.RDS',
        training_accuracy = 'outputs/singlem_abundtrim_kmer_optimal_rf/{study}_training_acc.csv',
        training_confusion = 'outputs/singlem_abundtrim_kmer_optimal_rf/{study}_training_confusion.pdf',
        validation_accuracy = 'outputs/singlem_abundtrim_kmer_optimal_rf/{study}_validation_acc.csv',
        validation_confusion = 'outputs/singlem_abundtrim_kmer_optimal_rf/{study}_validation_confusion.pdf'
    resources: 
        mem_mb=4000
    threads: 4
    params:
        threads = 4,
        validation_study = "{study}"
    conda: 'envs/tuneranger.yml'
    script: "scripts/tune_rf.R"


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

rule gather_gtdb_rep_to_shared_assemblies:
    input:  gather=expand("outputs/gather/{study}_vita_vars_gtdb_seed{seed}.csv", study = STUDY, seed = SEED) 
    output: gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv"
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/gather_gtdb_rep_to_shared_assemblies.R"


# TR TODO: SUMMARIZE TO SPECIES
# see: sandbox/test_gather_lineage_summarize/README.sh

#############################################
# Spacegraphcats Genome Queries
#############################################

# specifying make_sgc_conf as the target downloads the genomes of interest,
# circumventing a bug in genome grist that prevents using the target
# download gather genomes. The sgc conf file is a dummy file -- it will be
# written to outputs/sgc, but the conf file has the wrong catlas bases.
checkpoint download_shared_assemblies:
    input: 
        gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        conf = ""
    output: directory("genbank_genomes")
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run conf.yml --until make_sgc_conf
    '''

def aggregate_download_shared_assemblies:
    checkpoint_output = checkpoints.download_shared_assemblies.get(**wildcards).output[0]  
    file_names = expand("genbank_genomes/{genome}_genomic.fna.gz", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}_genomic.fna.gz")).genome)
    return file_names

rule make_sgc_conf_files:
    input:
        csv = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv"
        queries = "genbank_genomes/{acc}_genomic.fna.gz",
    output:
        conf = "outputs/sgc_conf/{library}_k31_r1.conf"
    run:
        query_list = "\n- ".join(input.queries)
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.library}
input_sequences:
- outputs/abundtrim/{wildcards.library}.abundtrim.fq.gz
ksize: 31
radius: 1
search:
- {query_list}
""", file=fp)

rule spacegraphcats_shared_assemblies:
    input: 
        query = "outputs/gather_matches_loso/{gather_genome}.gz", 
        conf = "inputs/sgc_conf/{library}_r1_conf.yml",
        reads = "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output:
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}.gz.cdbg_ids.reads.gz",
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}.gz.contigs.sig"
    params: outdir = "outputs/sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 64000
    threads: 1
    shell:'''
    spacegraphcats {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir}  
    '''

rule prokka_gather_match_genomes:
    output: 
        ffn = 'outputs/gather_shared_assemblies_prokka/{acc}.ffn',
        faa = 'outputs/gather_shared_assemblies_prokka/{acc}.faa'
    input: 'genbank_genomes/{acc}_genomic.fna.gz'
    conda: 'envs/prokka.yml'
    resources:
        mem_mb = 8000
    threads: 2
    params: 
        outdir = 'outputs/gather_shared_assemblies_prokka/',
        prefix = lambda wildcards: wildcards.acc,
        gzip = lambda wildcards: "outputs/gather_shared_assemblies/" + wildcards.acc + "_genome.fna"
    shell:'''
    gunzip {input}
    prokka {params.gzip} --outdir {params.outdir} --prefix {params.prefix} --metagenome --force --locustag {params.prefix} --cpus {threads} --centre X --compliant
    mv {params.prefix}.ffn {output.ffn}
    mv {params.prefix}.faa {output.faa}
    gzip {params.gzip}
    '''

# TR TODO: UPDATE ENV? 
# TR TODO: UPDATE INPUT FOR CHECKPOINT
#rule spacegraphcats_multifasta:
#    input:
#        queries = expand('outputs/gather_matches_loso_prokka/{acc}.ffn', gather_genome = GATHER_GENOMES),
#        conf = "inputs/sgc_conf/{library}_r1_multifasta_conf.yml",
#        reads = "outputs/abundtrim/{library}.abundtrim.fq.gz",
#        #sig = "outputs/vita_rf_seed/at_least_5_studies_vita_vars.sig"
#    output: "outputs/sgc_genome_queries/{library}_k31_r1_multifasta/query-results.csv"
#    params: 
#        outdir = "outputs/sgc_genome_queries",
#        #out = lambda wildcards: wildcards.library + "_k31_r1_multifasta/query-results.csv"
#    conda: "envs/spacegraphcats_multifasta.yml"
#    resources:
#        mem_mb = 32000
#    threads: 1
#    shell:'''
#    python -m spacegraphcats {input.conf} multifasta_query --nolock --outdir {params.outdir}
#    '''

##############################################
## Pangenome signature/variable importance
##############################################

# use default contig sigs output by sgc to start. 

rule calc_sig_nbhd_reads:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz"
    output: "outputs/nbhd_reads_sigs/{library}/{acc}.cdbg_ids.reads.sig"
    params: name = lambda wildcards: wildcards.library + "_" + wildcards.acc
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 2000
    threads: 1
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} --merge {params.name} {input}
    '''

rule nbhd_read_sig_to_csv:
    input: "outputs/nbhd_reads_sigs/{library}/{acc}.cdbg_ids.reads.sig"
    output: "outputs/nbhd_reads_sigs_csv/{library}/{acc}.cdbg_ids.reads.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 2000
    threads: 1
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

# TR TODO: UPDATE INPUT TO CHECKPOINT
#rule index_sig_nbhd_reads:
#    input: expand("outputs/nbhds_read_sigs/{library}/{acc}.cdbg_ids.reads.sig", library = LIBRARIES, gather_genome = GATHER_GENOMES)
#    output: "outputs/nbhd_reads_sigs_gather/nbhd_reads_sigs.sbt.json"
#    conda: "envs/sourmash.yml"
#    resources:
#        mem_mb = 2000
#    threads: 1
#    shell:'''
#    sourmash index -k 31 {output} --traverse-directory outputs/nbhd_read_sigs
#    '''

#rule gather_vita_vars_study_against_nbhd_read_sigs:
#    input:
#        sig="outputs/vita_rf_seed/at_least_5_studies_vita_vars.sig",
#        db = "outputs/nbhd_reads_sigs_gather/nbhd_read_sigs.sbt.json"
#    output: 
#        csv="outputs/nbhd_read_sigs_gather/at_least_5_studies_vita_vars_vs_nbhd_read_sigs_tbp0.csv",
#        un="outputs/nbhd_read_sigs_gather/at_least_5_studies_vita_vars_vs_nbhd_read_sigs_tbp0.un"
#    conda: 'envs/sourmash.yml'
#    resources:
#        mem_mb = 32000
#    threads: 1
#    shell:'''
#    sourmash gather -o {output.csv} --output-unassigned {output.un} --scaled 2000 --threshold-bp 0 -k 31 {input.sig} {input.db}
#    '''

rule merge_sgc_sigs_to_pangenome:
    input: expand("outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{{acc}}_genomic.fna.gz.contigs.sig", library = LIBRARIES)
    output: "outputs/sgc_pangenome_sigs/{acc}.sig"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 16000
    threads: 1
    shell:'''
    sourmash signature merge -o {output} -k 31 {input}
    '''

rule rename_sgc_pangenome_sigs:
    input: "outputs/sgc_pangenome_sigs/{acc}.sig"
    output: "outputs/sgc_pangenome_sigs/{acc}_renamed.sig"
    conda: "envs/sourmash.yml"
    resources:
        mem_mb = 1000
    threads: 1
    shell:'''
    sourmash sig rename -o {output} {input} {wildcards.gather_genome}_pangenome
    '''

# TR TODO: UPDATE INPUT TO CHECKPOINT
#rule compare_sgc_pangenome_sigs:
#    input: expand("outputs/sgc_pangenome_sigs/{acc}_renamed.sig", gather_genome = GATHER_GENOMES)
#    output: 
#        comp="outputs/sgc_pangenome_compare/pangenome_compare.comp",
#        csv="outputs/sgc_pangenome_compare/pangenome_compare.csv"
#    conda: "envs/sourmash.yml"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    sourmash compare -o {output.comp} --csv {output.csv} --ignore-abundance {input}
#    '''

# rule plot_sgc_pangenome_sigs:
#    input: "outputs/sgc_pangenome_compare/pangenome_compare.comp"
#    output: "outputs/sgc_pangenome_compare/pangenome_compare.comp.matrix.pdf"
#    conda:"envs/sourmash.yml"
#    resources:
#        mem_mb = 4000
#    threads: 1
#    shell:'''
#    sourmash plot --labels {input}
#    '''

# TR TODO: UPDATE INPUT TO CHECKPOINT
#rule index_sgc_pangenome_sigs:
#    input: expand("outputs/sgc_pangenome_sigs/{gather_genome}_renamed.sig", gather_genome = GATHER_GENOMES)
#    output: "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json"
#    conda: "envs/sourmash.yml"
#    resources:
#        mem_mb = 32000
#    threads: 1
#    shell:'''
#    sourmash index -k 31 {output} {input}
#    '''

#rule gather_vita_vars_study_against_sgc_pangenome_sigs:
#    input:
#        sig="outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.sig",
#        db = "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json"
#    output: 
#        csv="outputs/gather_sgc_pangenome/{study}_vita_vars_seed{seed}_pangenome.csv",
#        un="outputs/gather_sgc_pangenome/{study}_vita_vars_seed{seed}_pangenome.un"
#    conda: 'envs/sourmash.yml'
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    sourmash gather --threshold-bp 0 -o {output.csv} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
#    '''

# TR TODO: UPDATE INPUT TO GTDB DB
#rule gather_against_sgc_pangenome_sigs_plus_all_dbs:
#    input:
#        sig="outputs/vita_rf_seed/{study}_vita_vars_seed{seed}.sig",
#        db0 = "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json",
#        db1="/home/irber/sourmash_databases/outputs/sbt/genbank-bacteria-x1e6-k31.sbt.zip",
#        db2="/home/irber/sourmash_databases/outputs/sbt/genbank-viral-x1e6-k31.sbt.zip",
#        db3="/home/irber/sourmash_databases/outputs/sbt/genbank-archaea-x1e6-k31.sbt.zip",
#        db4="/home/irber/sourmash_databases/outputs/sbt/genbank-fungi-x1e6-k31.sbt.zip",
#        db5="/home/irber/sourmash_databases/outputs/sbt/genbank-protozoa-x1e6-k31.sbt.zip",
#    output: 
#        csv="outputs/sgc_pangenome_gather/{study}_vita_vars_all_seed{seed}.csv",
#        un="outputs/sgc_pangenome_gather/{study}_vita_vars_all_seed{seed}.un"
#    conda: 'envs/sourmash.yml'
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    sourmash gather -o {output.csv} --threshold-bp 0 --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db0} {input.db1} {input.db2} {input.db3} {input.db4} {input.db5}
#    '''

##############################################
## Query by multifasta eggnog gene annotation
##############################################

#rule combine_gather_match_genomes_prokka: 
#    input: expand("outputs/gather_matches_loso_prokka/{gather_genome}.faa", gather_genome = GATHER_GENOMES)
#    output: "outputs/gather_matches_loso_prokka/all_gather_genome_matches.faa"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    cat {input} > {output}
#    '''

#rule combine_multifasta_results:
#    input: expand("outputs/sgc_genome_queries/{library}_k31_r1_multifasta/query-results.csv", library = LIBRARIES)
#    output: "outputs/gather_matches_loso_multifasta/all-multifasta-query-results.csv"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    cat {input} > {output}
#    '''

#rule get_multifasta_query_gene_names:
#    input: "outputs/gather_matches_loso_multifasta/all-multifasta-query-results.csv"
#    output: names = "outputs/gather_matches_loso_multifasta/all-multifasta-query-results-names.txt"
#    conda: 'envs/tidy.yml'
#    resources:
#        mem_mb = 8000
#    threads: 1
#    script: "scripts/get_multifasta_names.R"

#rule grab_multifasta_query_genes:
#    output: "outputs/gather_matches_loso_multifasta/all-multifasta-query-results.faa"
#    input:
#        names = "outputs/gather_matches_loso_multifasta/all-multifasta-query-results-names.txt",
#        faa = 'outputs/gather_matches_loso_prokka/all_gather_genome_matches.faa'
#    conda: "envs/sourmash.yml"
#    resources:
#        mem_mb = 4000
#    threads: 1
#    shell: '''
#    scripts/extract-aaseq-matches.py {input.names} {input.faa} > {output}
#    '''

#rule eggnog_multifasta_query_genes:
#    input:
#        faa = "outputs/gather_matches_loso_multifasta/all-multifasta-query-results.faa",
#        db = "inputs/eggnog_db/eggnog.db"
#    output: "outputs/gather_matches_loso_multifasta/all-multifasta-query-results.emapper.annotations"
#    resources:
#        mem_mb = 64000
#    threads: 8
#    params: 
#        outdir = "outputs/gather_matches_loso_multifasta/",
#        dbdir = "inputs/eggnog_db",
#        out_prefix = "all-multifasta-query-results"
#    conda: 'envs/eggnog.yml'
#    shell:'''
#    emapper.py --cpu {threads} -i {input.faa} --output {params.out_prefix} --output_dir {params.outdir} -m diamond -d none --tax_scope auto --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0 --override --temp_dir tmp/ -d bact --data_dir {params.dbdir}
#    '''

###########################################################
## SingleM Ribosomal gene abundance validation: SGC NBHDS
###########################################################

#rule rename_sgc_nbhd_reads:
#    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz"
#    output: "outputs/sgc_genome_queries_renamed/{library}/{acc}_renamed.fastq.gz"
#    resources:
#        mem_mb = 4000
#    threads: 1
#    shell:'''
#    zcat {input} | awk '{{print (NR%4 == 1) ? "@1_" ++i : $0}}' | gzip -c > {output}
#    '''
#
#rule singlem_default_nbhd_reads:
#    input: "outputs/sgc_genome_queries_renamed/{library}/{acc}_renamed.fastq.gz"
#    output: "outputs/singlem_sgc_genome_queries/{library}/{acc}_otu_default.csv"
#    conda: "envs/singlem.yml"
#    resources:
#        mem_mb = 16000
#    threads: 2
#    shell: '''
#    singlem pipe --sequences {input} --otu_table {output} --output_extras --threads {threads} --filter_minimum_nucleotide 36 --min_orf_length 36 --filter_minimum_protein 12 || touch {output}
#    touch {output} # creates output file for runs with no seq matches
#    '''
#
#rule download_singlem_16s_pkg:
#    output: "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg.tar.xz"
#    resources:
#        mem_mb = 1000
#    threads: 1
#    shell:'''
#    wget -O {output} https://github.com/wwood/singlem_extra_packages/raw/master/release1/4.40.2013_08_greengenes_97_otus.with_euks.spkg.tar.xz
#    '''
#
#rule untar_singlem_16s_pkg:
#    input:"inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg.tar.xz"
#    output: "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg/CONTENTS.json"
#    resources:
#        mem_mb = 1000
#    threads: 1
#    params: outdir = "inputs/singlem"
#    shell:'''
#    tar xf {input} -C {params.outdir}
#    '''
#
#rule singlem_16s_nbhd_reads:
#    input: 
#        seq = "outputs/sgc_genome_queries_renamed/{library}/{acc}_renamed.fastq.gz",
#        pkg =  "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg/CONTENTS.json"
#    output: "outputs/singlem_sgc_genome_queries/{library}/{acc}_otu_16s.csv"
#    conda: "envs/singlem.yml"
#    resources:
#        mem_mb = 16000
#    threads: 2
#    params: 
#        pkg_dir = "inputs/singlem/4.40.2013_08_greengenes_97_otus.with_euks.spkg"
#    shell: '''
#    singlem pipe --sequences {input.seq} --singlem_packages {params.pkg_dir} --otu_table {output} --output_extras --threads {threads} --filter_minimum_nucleotide 36 --min_orf_length 36 --filter_minimum_protein 12 || touch {output}
#    touch {output}
#    '''

# TR TODO: UPDATE TO CHECKPOINT INPUT
#rule combine_singlem_default:
#    input:
#        default = expand("outputs/singlem_sgc_genome_queries/{library}/{acc}_otu_default.csv", library = LIBRARIES, gather_genome = GATHER_GENOMES),
#    output: res = "outputs/singlem_sgc_genome_queries/combined_default.tsv"
#    conda: "envs/tidy.yml"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    script: "scripts/parse_singlem_default.R"

# TR TODO: UPDATE TO CHECKPOINT INPUT
#rule combine_singlem_16s:
#    input:
#        s16 = expand("outputs/singlem_sgc_genome_queries/{library}/{acc}_otu_16s.csv", library = LIBRARIES, gather_genome = GATHER_GENOMES),
#    output: res = "outputs/singlem_sgc_genome_queries/combined_16s.tsv"
#    conda: "envs/tidy.yml"
#    resources:
#        mem_mb = 8000
#    threads: 1
#    script: "scripts/parse_singlem_16s.R"

# rule combine_singlem:
#     input: 
#         s16 = "outputs/singlem_sgc_genome_queries/combined_16s.tsv",
#         default = "outputs/singlem_sgc_genome_queries/combined_default.tsv"
#     output: res = "outputs/singlem_sgc_genome_queries/combined.tsv"
#     conda: "envs/tidy.yml"
#     resources:
#         mem_mb = 8000
#     threads: 1
#     script: "scripts/parse_singlem.R"
# 
# rule singlem_to_counts:
#     input:  res = "outputs/singlem_sgc_genome_queries/combined.tsv"
#     output: counts = "outputs/singlem_sgc_genome_queries/singlem_counts.tsv"
#     conda: "envs/tidy.yml"
#     resources:
#         mem_mb = 8000
#     threads: 1
#     script: "scripts/make_singlem_counts.R"
# 
# rule singlem_install_pomona:
#     input: "outputs/singlem_sgc_genome_queries/combined.tsv"
#     output:
#         pomona = "outputs/singlem_vita_rf/pomona_install.txt"
#     resources:
#         mem_mb = 1000
#     threads: 1
#     conda: 'envs/rf.yml'
#     script: "scripts/install_pomona.R"
# 
# rule singlem_var_sel_rf:
#     input:
#         info = "inputs/working_metadata.tsv", 
#         counts = "outputs/singlem_sgc_genome_queries/singlem_counts.tsv",
#         pomona = "outputs/singlem_vita_rf/pomona_install.txt"
#     output:
#         vita_rf = "outputs/singlem_vita_rf/{study}_vita_rf.RDS",
#         vita_vars = "outputs/singlem_vita_rf/{study}_vita_vars.txt",
#         ibd_filt = "outputs/singlem_vita_rf/{study}_ibd_filt.csv"
#     resources:
#         mem_mb = 64000
#     threads: 8
#     params: 
#         threads = 8,
#         validation_study = "{study}"
#     conda: 'envs/rf.yml'
#     script: "scripts/singlem_vita_rf.R"
# 
# 
# rule singlem_loo_validation:
#     input: 
#         ibd_filt = 'outputs/singlem_vita_rf/{study}_ibd_filt.csv',
#         info = 'inputs/working_metadata.tsv',
#         eval_model = 'scripts/function_evaluate_model.R',
#         ggconfusion = 'scripts/ggplotConfusionMatrix.R'
#     output: 
#         recommended_pars = 'outputs/singlem_optimal_rf/{study}_rec_pars.tsv',
#         optimal_rf = 'outputs/singlem_optimal_rf/{study}_optimal_rf.RDS',
#         training_accuracy = 'outputs/singlem_optimal_rf/{study}_training_acc.csv',
#         training_confusion = 'outputs/singlem_optimal_rf/{study}_training_confusion.pdf',
#         validation_accuracy = 'outputs/singlem_optimal_rf/{study}_validation_acc.csv',
#         validation_confusion = 'outputs/singlem_optimal_rf/{study}_validation_confusion.pdf'
#     resources:
#         mem_mb = 8000
#     threads: 8
#     params:
#         threads = 8,
#         validation_study = "{study}"
#     conda: 'envs/tuneranger.yml'
#     script: "scripts/singlem_tune_rf.R"
# 
# rule extract_singlem_read_names_default_sgc:
#     input: singlem = "outputs/singlem_sgc_genome_queries/{library}/{gather_genome}_otu_default.csv",
#     output: reads = "outputs/singlem_sgc_genome_queries_reads/{library}/{gather_genome}_otu_default_names.txt" 
#     conda: "envs/tidy.yml"
#     resources:
#         mem_mb = 8000
#     threads: 1
#     script: "scripts/extract_singlem_read_names.R"
# 
# rule extract_singlem_reads_default_sgc:
#     input: 
#         names = "outputs/singlem_sgc_genome_queries_reads/{library}/{gather_genome}_otu_default_names.txt",
#         fq = "outputs/sgc_genome_queries_renamed/{library}/{gather_genome}_renamed.fastq.gz",
#     output: "outputs/singlem_sgc_genome_queries_reads/{library}/{gather_genome}_otu_default.fq",
#     conda: "envs/sourmash.yml"
#     resources:
#         mem_mb = 8000
#     threads: 1
#     shell:'''
#     scripts/extract-aaseq-matches.py {input.names} {input.fq} > {output}
#     '''
# 
# rule extract_singlem_read_names_16s_sgc:
#     input: singlem = "outputs/singlem_sgc_genome_queries/{library}/{gather_genome}_otu_16s.csv",
#     output: reads = "outputs/singlem_sgc_genome_queries_reads/{library}/{gather_genome}_otu_16s_names.txt" 
#     conda: "envs/tidy.yml"
#     resources:
#         mem_mb = 8000
#     threads: 1
#     script: "scripts/extract_singlem_read_names.R"
# 
# rule extract_singlem_reads_16s_sgc:
#     input: 
#         names = "outputs/singlem_sgc_genome_queries_reads/{library}/{gather_genome}_otu_16s_names.txt",
#         fq = "outputs/sgc_genome_queries_renamed/{library}/{gather_genome}_renamed.fastq.gz",
#     output: "outputs/singlem_sgc_genome_queries_reads/{library}/{gather_genome}_otu_16s.fq",
#     resources:
#         mem_mb = 8000
#     threads: 1
#     conda: "envs/sourmash.yml"
#     shell:'''
#     scripts/extract-aaseq-matches.py {input.names} {input.fq} > {output}
#     '''
# 
# rule combine_singlem_reads_per_lib_sgc:
#     input:
#         expand("outputs/singlem_sgc_genome_queries_reads/{{library}}/{gather_genome}_otu_default.fq", gather_genome = GATHER_GENOMES),
#         expand("outputs/singlem_sgc_genome_queries_reads/{{library}}/{gather_genome}_otu_16s.fq", gather_genome = GATHER_GENOMES)
#     output: "outputs/singlem_sgc_genome_queries_reads/{library}_singlem_reads.fq"
#     resources:
#         mem_mb = 8000
#     threads: 1
#     shell:'''
#     cat {input} > {output}
#     '''
# 
# ## run random forests on signatures of singlem output ###########################
#   
# rule compute_signatures_singlem_sgc:
#     input: "outputs/singlem_sgc_genome_queries_reads/{library}_singlem_reads.fq"
#     output: "outputs/singlem_sgc_genome_queries_sigs/{library}_singlem_reads.sig"
#     conda: "envs/sourmash.yml"
#     resources:
#         mem_mb = 1000
#     threads: 1
#     shell:'''
#     sourmash compute -k 31 --scaled 2000 -o {output} --track-abundance {input} || touch {output} 
#     '''
# 
# rule convert_greater_than_1_signatures_to_csv_singlem_sgc:
#     input: "outputs/singlem_sgc_genome_queries_sigs/{library}_singlem_reads.sig"
#     output: "outputs/singlem_sgc_genome_queries_kmer_csv/{library}_singlem_reads.csv"
#     conda: 'envs/sourmash.yml'
#     resources: 
#         mem_mb=1000
#     threads: 1
#     shell:'''
#     python scripts/sig_to_csv.py {input} {output} || touch {output}
#     '''
# 
# rule make_hash_abund_table_long_normalized_singlem_sgc:
#     input: 
#         expand("outputs/singlem_sgc_genome_queries_kmer_csv/{library}_singlem_reads.csv", library = LIBRARIES)
#     output: csv = "outputs/singlem_sgc_genome_queries_kmer_hash_tables/normalized_abund_hashes_long.csv"
#     conda: 'envs/r.yml'
#     resources: 
#         mem_mb=16000
#     threads: 1
#     script: "scripts/normalized_hash_abund_long_singlem.R"
# 
# rule make_hash_abund_table_wide_singlem_sgc:
#     input: "outputs/singlem_sgc_genome_queries_kmer_hash_tables/normalized_abund_hashes_long.csv"
#     output: "outputs/singlem_sgc_genome_queries_kmer_hash_tables/normalized_abund_hashes_wide.feather"
#     resources: 
#         mem_mb=32000
#     threads: 1
#     run:
#         import pandas as pd
#         import feather
#         
#         ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
#         ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
#         ibd_wide = ibd_wide.fillna(0)
#         ibd_wide['sample'] = ibd_wide.index
#         ibd_wide = ibd_wide.reset_index(drop=True)
#         ibd_wide.columns = ibd_wide.columns.astype(str)
#         ibd_wide.to_feather(str(output)) 
# 
# rule singlem_kmer_install_pomona_sgc:
#     input: "outputs/singlem_sgc_genome_queries_kmer_hash_tables/normalized_abund_hashes_wide.feather"
#     output:
#         pomona = "outputs/singlem_sgc_genome_queries_kmer_vita_rf/pomona_install.txt"
#     conda: 'envs/rf.yml'
#     resources: 
#         mem_mb=1000
#     threads: 1
#     script: "scripts/install_pomona.R"
# 
# rule singlem_kmer_vita_var_sel_rf_sgc:
#     input:
#         info = "inputs/working_metadata.tsv", 
#         feather = "outputs/singlem_sgc_genome_queries_kmer_hash_tables/normalized_abund_hashes_wide.feather",
#         pomona = "outputs/singlem_sgc_genome_queries_kmer_vita_rf/pomona_install.txt"
#     output:
#         vita_rf = "outputs/singlem_sgc_genome_queries_kmer_vita_rf/{study}_vita_rf.RDS",
#         vita_vars = "outputs/singlem_sgc_genome_queries_kmer_vita_rf/{study}_vita_vars.txt",
#         ibd_filt = "outputs/singlem_sgc_genome_queries_kmer_vita_rf/{study}_ibd_filt.csv"
#     resources: 
#         mem_mb=16000
#     threads: 4
#     params: 
#         threads = 4,
#         validation_study = "{study}"
#     conda: 'envs/rf.yml'
#     script: "scripts/singlem_kmer_vita_rf.R"
# 
# rule singlem_kmer_loo_validation_sgc:
#     input: 
#         ibd_filt = 'outputs/singlem_sgc_genome_queries_kmer_vita_rf/{study}_ibd_filt.csv',
#         info = 'inputs/working_metadata.tsv',
#         eval_model = 'scripts/function_evaluate_model.R',
#         ggconfusion = 'scripts/ggplotConfusionMatrix.R'
#     output: 
#         recommended_pars = 'outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_rec_pars.tsv',
#         optimal_rf = 'outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_optimal_rf.RDS',
#         training_accuracy = 'outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_training_acc.csv',
#         training_confusion = 'outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_training_confusion.pdf',
#         validation_accuracy = 'outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_validation_acc.csv',
#         validation_confusion = 'outputs/singlem_sgc_genome_queries_kmer_optimal_rf/{study}_validation_confusion.pdf'
#     resources: 
#         mem_mb=4000
#     threads: 4
#     params:
#         threads = 4,
#         validation_study = "{study}"
#     conda: 'envs/tuneranger.yml'
#     script: "scripts/tune_rf.R"

##############################################
## Pangenome differential abundance analysis
##############################################

#rule diginorm_nbhd_reads:
#    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.gz"
#    output: "outputs/nbhd_reads_diginorm/{library}/{gather_genome}.cdgb_ids.reads.diginorm.gz"
#    conda: "envs/env.yml"
#    resources:
#        mem_mb = 16000
#    threads: 1
#    shell:'''
#    normalize-by-median.py -k 20 -C 20 -M 16e9 --gzip -o {output} {input}
#    '''
#
## TR TODO: UPDATE MEGAHIT TO USE PAIRED-END READS
#rule megahit:
#    input: "outputs/nbhd_reads_diginorm/{library}/{gather_genome}.cdgb_ids.reads.diginorm.gz"
#    output: "outputs/nbhd_reads_diginorm_megahit/{library}/{gather_genome}.contigs.fa"
#    conda: 'envs/assembly.yml'
#    resources: 
#        mem_mb = 8000
#    threads: 2
#    shell:'''
#    megahit -r {input} -t {threads} --min-contig-len 500 \
#        --out-dir {wildcards.library}_{wildcards.gather_genome}_megahit \
#        --out-prefix {wildcards.library}_{wildcards.gather_genome}
#    mv  {wildcards.library}_{wildcards.gather_genome}_megahit/{wildcards.library}_{wildcards.gather_genome}.contigs.fa {output}
#    rm -rf {wildcards.library}_{wildcards.gather_genome}_megahit
#    '''
#
## calculate percent of reads from each nbhd that map
#
#rule index:
#    input: genome = "outputs/nbhd_reads_diginorm_megahit/{library}/{gather_genome}.contigs.fa"
#    output:  "outputs/nbhd_reads_diginorm_megahit/{library}/{gather_genome}.contigs.fa.bwt"
#    conda: "envs/assembly.yml"
#    resources: mem_mb = 8000
#    threads: 1
#    shell:'''
#    bwa index {input}
#    '''
#
## TR TODO: UPDATE TO PAIRED END
#rule bwa:
#    input: 
#        indx =  "outputs/nbhd_reads_diginorm_megahit/{library}/{gather_genome}.contigs.fa.bwt",
#        genome = "outputs/nbhd_reads_diginorm_megahit/{library}/{gather_genome}.contigs.fa",
#        reads =  "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.gz"
#    output: "outputs/nbhd_reads_diginorm_megahit_bwa/{library}/{gather_genome}.bam"
#    conda: "envs/assembly.yml"
#    resources: mem_mb = 3000
#    threads: 2
#    shell:'''
#    bwa mem -t {threads} {input.genome} {input.reads} | samtools sort -o {output} - || touch {output}
#    ''' 
#
#rule samtools_flagstat:
#    input: "outputs/nbhd_reads_diginorm_megahit_bwa/{library}/{gather_genome}.bam"
#    output: "outputs/nbhd_reads_diginorm_megahit_bwa/{library}/{gather_genome}.flagstat"
#    conda: "envs/assembly.yml"
#    resources: mem_mb = 2000
#    threads: 1
#    shell:'''
#    samtools flagstat {input} > {output} - || touch {output}
#    '''
#
## check alignment to single assembly first; 
## then use prokka to predict ORFs, 
## cdhit ORF sequences,
## and align all reads to pangenome sequence. 
#rule prokka_megahit:
#    output: 
#        ffn = 'outputs/nbhd_reads_diginorm_megahit_prokka/{library}/{gather_genome}.ffn',
#        faa = 'outputs/nbhd_reads_diginorm_megahit_prokka/{library}/{gather_genome}.faa'
#    input: 'outputs/nbhd_reads_diginorm_megahit/{library}/{gather_genome}.contigs.fa'
#    conda: 'envs/assembly.yml'
#    resources: mem_mb = 4000
#    threads: 2
#    params: 
#        output_folder = lambda wildcards: 'prokka/' + wildcards.library
#    shell:'''
#    prokka {input} --outdir {params.output_folder} --prefix {wildcards.gather_genome} --metagenome --force --locustag {wildcards.library} --cpus {threads} || touch {output.ffn}
#    touch {output.faa}
#    '''
#
#rule cat_prokka:
#    input: expand('outputs/nbhd_reads_diginorm_megahit_prokka/{library}/{{gather_genome}}.ffn', library = LIBRARIES)
#    output: 'outputs/nbhd_reads_diginorm_megahit_cat_prokka/{gather_genome}.ffn'
#    resources: mem_mb = 2000
#    threads: 1
#    shell:'''
#    cat {input} > {output}
#    '''
#
#rule cdhit:
#    input: 'outputs/nbhd_reads_diginorm_megahit_cat_prokka/{gather_genome}.ffn'
#    output: 'outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit/{gather_genome}.cdhit.ffn'
#    conda: 'envs/assembly.yml'
#    resources: mem_mb = 16000
#    threads: 2
#    shell:'''
#    cd-hit-est -i {input} -o {output} -c 0.9 -n 8 -d 200 -M 15500 -T {threads}
#    '''
#
#rule index_cdhit:
#    input: genome = "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit/{gather_genome}.cdhit.ffn"
#    output:  "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit/{gather_genome}.cdhit.ffn.bwt"
#    conda: "envs/assembly.yml"
#    resources: mem_mb = 2000
#    threads: 1
#    shell:'''
#    bwa index {input}
#    '''
#
## TR TODO: UPDATE TO PAIRED-END?
#rule bwa_cdhit:
#    input: 
#        indx =  "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit/{gather_genome}.cdhit.ffn.bwt",
#        genome = "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit/{gather_genome}.cdhit.ffn",
#        reads =  "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.gz"
#    output: "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit_bwa/{library}/{gather_genome}.bam"
#    conda: "envs/assembly.yml"
#    threads: 2
#    resources: mem_mb = 2000
#    shell:'''
#    bwa mem -t {threads} {input.genome} {input.reads} | samtools sort -o {output} -
#    ''' 
#
#rule flagstat_cdhit:
#    input: "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit_bwa/{library}/{gather_genome}.bam"
#    output: "outputs/nbhd_reads_diginorm_megahit_cat_prokka_cdhit_bwa/{library}/{gather_genome}.flagstat"
#    conda: "envs/assembly.yml"
#    shell:'''
#    samtools flagstat {input} > {output}
#    '''

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
