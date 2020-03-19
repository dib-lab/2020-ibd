import pandas as pd
import feather
from sourmash import signature
import glob
import os
from collections import Counter

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
SAMPLES = m.sort_values(by='read_count')['run_accession']
LIBRARIES = m['library_name'].unique().tolist()
STUDY = m['study_accession'].unique().tolist()

rule all:
    input:
        # sourmash compare outputs:
        "outputs/comp/all_filt_permanova_cosine.csv",
        "outputs/comp/all_filt_permanova_jaccard.csv",
        "outputs/comp/study_plt_all_filt_jaccard.pdf",
        "outputs/comp/diagnosis_plt_all_filt_jaccard.pdf",
        "outputs/comp/study_plt_all_filt_cosine.pdf",
        "outputs/comp/diagnosis_plt_all_filt_cosine.pdf",
        # variable selection outputs:
        "outputs/filt_sig_hashes/count_total_hashes.txt",
        expand("outputs/vita_rf/{study}_vita_rf.RDS", study = STUDY),
        expand("outputs/vita_rf/{study}_vita_vars.txt", study = STUDY),
        expand("outputs/vita_rf/{study}_ibd_filt.csv", study = STUDY),
        # optimal RF outputs:
        expand('outputs/optimal_rf/{study}_optimal_rf.RDS', study = STUDY),
        # variable characterization outputs
        expand("outputs/gather/{study}_vita_vars_refseq.csv", study = STUDY),
        expand("outputs/gather/{study}_vita_vars_genbank.csv", study = STUDY),
        expand("outputs/gather/{study}_vita_vars_all.csv", study = STUDY)
        #"outputs/gtdbtk/gtdbtk.bac120.summary.tsv",
        #expand("outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genomes}.fna.contigs.sig", 
        #       library = LIBRARIES, gather_genomes = GATHER_GENOMES)
        #"aggregated_checkpoints/aggregate_spacegraphcats_gather_matches.txt",
        #"aggregated_checkpoints/aggregate_spacegraphcats_gather_matches_plass.txt",
        #"outputs/hash_tables/all_unnormalized_abund_hashes_wide.feather",

########################################
## PREPROCESSING
########################################

rule download_fastq_files_R1:
    output: 
        r1="inputs/raw/{sample}_1.fastq.gz",
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastq_1 = row['fastq_ftp_1'].values
        fastq_1 = fastq_1[0]
        shell("wget -O {output.r1} {fastq_1}")


rule download_fastq_files_R2:
    output:
        r2="inputs/raw/{sample}_2.fastq.gz"
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastq_2 = row['fastq_ftp_2'].values
        fastq_2 = fastq_2[0]
        shell("wget -O {output.r2} {fastq_2}")


rule cat_libraries_R1:
    input: expand("inputs/raw/{sample}_1.fastq.gz", sample = SAMPLES)
    output: expand("inputs/cat/{library}_1.fastq.gz", library = LIBRARIES)
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
    conda: 'env.yml'
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
    conda: 'env2.yml'
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
    conda: 'env2.yml'
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
    conda: 'env.yml'
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{library}_R1.nohost.fq.gz',
        'outputs/bbduk/{library}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    conda: 'env.yml'
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule compute_signatures:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: "outputs/sigs/{library}.sig"
    conda: 'env.yml'
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} {input}
    '''

########################################
## Filtering and formatting signatures
########################################

rule get_greater_than_1_filt_sigs:
    input: expand("outputs/sigs/{library}.sig", library = LIBRARIES) 
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    run:
        # Determine the number of hashes, the number of unique hashes, and the number of
        # hashes that occur once across 954 IBD/control gut metagenomes (excludes the 
        # iHMP). Calculated for a scaled of 2k. 9 million hashes is the current 
        # approximate upper limit with which to build a sample vs hash abundance table 
        # using my current methods.

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
    input: expand("outputs/sigs/{library}.sig", library = LIBRARIES)
    output: "outputs/filt_sig_hashes/count_total_hashes.txt"
    run:
        files = input

        all_mins = []
        for file in files:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[1]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes
                all_mins += mins

        with open(str(output), "w") as f:
            print(len(all_mins), file=f)


rule convert_greater_than_1_hashes_to_sig:
    input: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name greater_than_one_count_hashes --filename {input} {input}
    '''

rule filter_signatures_to_greater_than_1_hashes:
    input:
        filt_sig = "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig",
        sigs = "outputs/sigs/{library}.sig"
    output: "outputs/filt_sigs/{library}_filt.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_filtered_sigs:
    input: "outputs/filt_sigs/{library}_filt.sig"
    output: "outputs/filt_sigs_named/{library}_filt_named.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.library}_filt
    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/filt_sigs_named/{library}_filt_named.sig"
    output: "outputs/filt_sigs_named_csv/{library}_filt_named.csv"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/filt_sigs_named_csv/{library}_filt_named.csv", library = LIBRARIES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'r.yml'
    script: "scripts/normalized_hash_abund_long.R"

rule make_hash_abund_table_wide:
    input: "outputs/hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
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
        pomona = "outputs/vita_rf/pomona_install.txt"
    conda: 'rf.yml'
    script: "scripts/install_pomona.R"

rule vita_var_sel_rf:
    input:
        info = "inputs/working_metadata.tsv", 
        feather = "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/vita_rf/{study}_vita_rf.RDS",
        vita_vars = "outputs/vita_rf/{study}_vita_vars.txt",
        ibd_filt = "outputs/vita_rf/{study}_ibd_filt.csv"
    params: 
        threads = 32,
        validation_study = "{study}"
    conda: 'rf.yml'
    script: "scripts/vita_rf.R"

rule loo_validation:
    input: 
        ibd_filt = 'outputs/vita_rf/{study}_ibd_filt.csv',
        info = 'inputs/working_metadata.tsv',
        eval_model = 'scripts/function_evaluate_model.R',
        ggconfusion = 'scripts/ggplotConfusionMatrix.R'
    output: 
        recommended_pars = 'outputs/optimal_rf/{study}_rec_pars.tsv',
        optimal_rf = 'outputs/optimal_rf/{study}_optimal_rf.RDS',
        training_accuracy = 'outputs/optimal_rf/{study}_training_acc.csv',
        training_confusion = 'outputs/optimal_rf/{study}_training_confusion.pdf',
        validation_accuracy = 'outputs/optimal_rf/{study}_validation_acc.csv',
        validation_confusion = 'outputs/optimal_rf/{study}_validation_confusion.pdf'
    params:
        threads = 20,
        validation_study = "{study}"
    conda: 'tuneranger.yml'
    script: "scripts/tune_rf.R"


############################################
## Predictive hash characterization - gather
############################################

rule convert_vita_vars_to_sig:
    input: "outputs/vita_rf/{study}_vita_vars.txt"
    output: "outputs/vita_rf/{study}_vita_vars.sig"
    conda: "sourmash.yml"
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name vita_vars --filename {input} {input}
    '''

rule download_gather_almeida:
    output: "inputs/gather_databases/almeida-mags-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/5jyzr/download
    '''

rule untar_almeida:
    output: "inputs/gather_databases/almeida-mags-k31.sbt.json"
    input: "inputs/gather_databases/almeida-mags-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_pasolli:
    output: "inputs/gather_databases/pasolli-mags-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/3vebw/download
    '''

rule untar_pasolli:
    output: "inputs/gather_databases/pasolli-mags-k31.sbt.json"
    input: "inputs/gather_databases/pasolli-mags-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_nayfach:
    output: "inputs/gather_databases/nayfach-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/y3vwb/download
    '''

rule untar_nayfach:
    output: "inputs/gather_databases/nayfach-k31.sbt.json"
    input: "inputs/gather_databases/nayfach-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.tar.gz"
    shell:'''
    wget -O {output} https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k31.tar.gz
    '''

rule untar_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.sbt.json"
    input:  "inputs/gather_databases/genbank-d2-k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_refseq:
    output: "inputs/gather_databases/refseq-d2-k31.tar.gz"
    shell:'''
    wget -O {output} https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k31.tar.gz
    '''

rule untar_refseq:
    output: "inputs/gather_databases/refseq-d2-k31.sbt.json"
    input:  "inputs/gather_databases/refseq-d2-k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule gather_vita_vars_all:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json"
    output: 
        csv="outputs/gather/{study}_vita_vars_all.csv",
        matches="outputs/gather/{study}_vita_vars_all.matches",
        un="outputs/gather/{study}_vita_vars_all.un"
    conda: 'sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db1} {input.db4} {input.db3} {input.db2}
    '''

rule gather_vita_vars_genbank:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db="inputs/gather_databases/genbank-d2-k31.sbt.json",
    output: 
        csv="outputs/gather/{study}_vita_vars_genbank.csv",
        matches="outputs/gather/{study}_vita_vars_genbank.matches",
        un="outputs/gather/{study}_vita_vars_genbank.un"
    conda: 'sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

rule gather_vita_vars_refseq:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db="inputs/gather_databases/refseq-d2-k31.sbt.json",
    output: 
        csv="outputs/gather/{study}_vita_vars_refseq.csv",
        matches="outputs/gather/{study}_vita_vars_refseq.matches",
        un="outputs/gather/{study}_vita_vars_refseq.un"
    conda: 'sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''


# THESE GENOMES WILL NEED TO BE UPDATED WITH THE NEW OUTPUT
# FROM VARIABLE SELECTION
#rule download_gather_match_genomes:
#    output: "outputs/gather/gather_genomes.tar.gz"
#    shell:'''
#    wget -O {output} https://osf.io/ungza/download
#    '''

#rule untar_gather_match_genomes:
#    output:  directory("outputs/gather_genomes/")
#    input:"outputs/gather/gather_genomes.tar.gz"
#    params: outdir = "outputs/gather_genomes"
#    shell:'''
#    mkdir -p {params.outdir}
#    tar xf {input} -C {params.outdir}
#    '''

# THIS RULE SHOULD BE REPLACED WITH SOURMASH LCA CLASSIFY
#rule gtdbtk_gather_matches:
#    """
#    this rule require the gtdbtk databases. The tool finds the database by 
#    using a path specified in a file in the environment. I predownloaded the 
#    databases and placed them in the required location.
#    The path is in this file:
#    .snakemake/conda/9de8946b/etc/conda/activate.d/gtdbtk.sh
#    """
#    input: directory("outputs/gather_genomes/")
#    #input: aggregate_decompress_gather_matches
#    output: "outputs/gtdbtk/gtdbtk.bac120.summary.tsv"
#    params:  outdir = "outputs/gtdbtk"
#    conda: "gtdbtk.yml"
#    shell:'''
#    gtdbtk classify_wf --genome_dir {input} --out_dir {params.outdir} --cpus 8 
#    '''

#############################################
# Spacegraphcats Genome Queries
#############################################

checkpoint spacegraphcats_gather_matches:
    input: 
        query = directory("outputs/gather_genomes/"),
        conf = expand("inputs/sgc_conf/{library}_r1_conf.yml", library = LIBRARIES),
        reads = expand("outputs/abundtrim/{library}.abundtrim.fq.gz", library = LIBRARIES)
    output: 
        directory(expand("outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/", library = LIBRARIES))
        #"outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.cdbg_ids.reads.fa.gz",
    params: outdir = "outputs/sgc_genome_queries"
    conda: "spacegraphcats.yml"
    shell:'''
    python -m spacegraphcats {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir}  
    '''

rule calc_sig_nbhd_reads:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.fna.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_read_sigs/{library}/{gather_genome}.cdbg_ids.reads.sig"
    conda: "sourmash.yml"
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} --merge {wildcards.library}_{wildcards.gather_genome} {input}
    '''

def aggregate_spacegraphcats_gather_matches(wildcards):
    # checkpoint_output produces the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.spacegraphcats_gather_matches.get(**wildcards).output[0]    
    file_names = expand("outputs/nbhd_read_sigs/{library}/{gather_genome}.cdbg_ids.reads.sig",
                        library = LIBRARIES, 
                        gather_genome = glob_wildcards(os.path.join(checkpoint_output, "{gather_genome}.fna.cdbg_ids.reads.fa.gz")).gather_genome)
    # file_names will return all 129 queries.
    # because this takes a long time, we will subset the file names returned
    # to the 4 nbhds that account for the largest number of predictive hashes.
    bacteroides = [f for f in file_names if "SRS476121_69" in f]
    faecalibacterium =  [f for f in file_names if "SRS147022_17" in f]
    rtorques =  [f for f in file_names if "GCA_001406235.1_14207_7_41_genomic" in f]
    fplautii =  [f for f in file_names if "GCA_001405435.1_14207_7_29_genomic" in f]
    select_file_names = bacteroides + faecalibacterium + rtorques + fplautii
    return select_file_names


rule aggregate_signatures:
    input: aggregate_spacegraphcats_gather_matches
    output: "aggregated_checkpoints/aggregate_spacegraphcats_gather_matches.txt"
    shell:'''
    touch {output}
    '''

rule plass_nbhd_reads:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.fna.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_read_plass/{library}/{gather_genome}.cdbg_ids.reads.plass.faa"
    conda: "plass.yml"
    shell:'''
    plass assemble {input} {output} tmp
    '''

rule cdhit_plass:
    input: "outputs/nbhd_read_plass/{library}/{gather_genome}.cdbg_ids.reads.plass.faa"
    output: "outputs/nbhd_read_cdhit/{library}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "plass.yml"
    shell:'''
    cd-hit -i {input} -o {output} -c 1
    '''

def aggregate_spacegraphcats_gather_matches_plass(wildcards):
    # checkpoint_output produces the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.spacegraphcats_gather_matches.get(**wildcards).output[0]    
    file_names = expand("outputs/nbhd_read_cdhit/{library}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa",
                        library = LIBRARIES, 
                        gather_genome = glob_wildcards(os.path.join(checkpoint_output, "{gather_genome}.fna.cdbg_ids.reads.fa.gz")).gather_genome)
    bacteroides = [f for f in file_names if "SRS476121_69" in f]
    faecalibacterium =  [f for f in file_names if "SRS147022_17" in f]
    rtorques =  [f for f in file_names if "GCA_001406235.1_14207_7_41_genomic" in f]
    fplautii =  [f for f in file_names if "GCA_001405435.1_14207_7_29_genomic" in f]
    select_file_names = bacteroides + faecalibacterium + rtorques + fplautii
    return select_file_names

rule aggregate_spacegraphcats_gather_matches_plass:
    input: aggregate_spacegraphcats_gather_matches_plass
    output: "aggregated_checkpoints/aggregate_spacegraphcats_gather_matches_plass.txt"
    shell:'''
    touch {output}
    '''

rule paladin_index_plass:
    input: "outputs/nbhd_read_cdhit/{library}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    output: "outputs/nbhd_read_cdhit/{library}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa.bwt"
    conda: "plass.yml"
    shell: '''
    paladin index -r3 {input}
    '''

rule paladin_align_plass:
    input:
        indx="outputs/nbhd_read_cdhit/{library}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa.bwt",
        reads="outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.fna.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_read_paladin/{library}/{gather_genome}.sam"
    params: indx = "outputs/nbhd_read_cdhit/{library}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "plass.yml"
    shell:'''
    paladin align -f 125 -t 2 {params.indx} {input.reads} > {output}
    '''

########################################
## PCoA
########################################

rule compare_signatures_cosine:
    input: 
        expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES),
    output: "outputs/comp/all_filt_comp_cosine.csv"
    conda: "sourmash.yml"
    shell:'''
    sourmash compare -k 31 -p 8 --csv {output} {input}
    '''

rule compare_signatures_jaccard:
    input: 
        expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES),
    output: "outputs/comp/all_filt_comp_jaccard.csv"
    conda: "sourmash.yml"
    shell:'''
    sourmash compare --ignore-abundance -k 31 -p 8 --csv {output} {input}
    '''

rule permanova_jaccard:
    input: 
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        perm = "outputs/comp/all_filt_permanova_jaccard.csv"
    conda: "vegan.yml"
    script: "scripts/run_permanova.R"

rule permanova_cosine:
    input: 
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        perm = "outputs/comp/all_filt_permanova_cosine.csv"
    conda: "vegan.yml"
    script: "scripts/run_permanova.R"

rule plot_comp_jaccard:
    input:
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        study = "outputs/comp/study_plt_all_filt_jaccard.pdf",
        diagnosis = "outputs/comp/diagnosis_plt_all_filt_jaccard.pdf"
    conda: "ggplot.yml"
    script: "scripts/plot_comp.R"

rule plot_comp_cosine:
    input:
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        study = "outputs/comp/study_plt_all_filt_cosine.pdf",
        diagnosis = "outputs/comp/diagnosis_plt_all_filt_cosine.pdf"
    conda: "ggplot.yml"
    script: "scripts/plot_comp.R"

########################################
## Differential abundance
########################################

rule hash_table_long_unnormalized:
    """
    Unlike the hashtable that is input into the random forest analysis, this
    hash table is not normalized by number of hashes in the filtered signature. 
    Differential expression software that we will be using to calculate 
    differential abundance expects unnormalized counts.
    """
    input: expand("outputs/filt_sigs_named_csv/{library}_filt_named.csv", library = LIBRARIES)
    output: csv = "outputs/hash_tables/unnormalized_abund_hashes_long.csv"
    conda: 'r.yml'
    script: "scripts/all_unnormalized_hash_abund_long.R"

rule hash_table_wide_unnormalized:
    input: "outputs/hash_tables/unnormalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/unnormalized_abund_hashes_wide.feather"
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

#rule differential_abundance_all:
#    input: "outputs/hash_tables/all_unnormalized_abund_hashes_wide.feather"

