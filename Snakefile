import pandas as pd
#import feather
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

# this variable is output after random forests. 
# Specifying this variable can be avoided by using checkpoints,
# but this makes the DAG take forever to solve. 
# this step requires user input anyway to download the matching
# random forest genomes, so specifying this variable manually is a 
# compromise. 
GATHER_GENOMES = ["ERS235530_10.fna", "ERS235531_43.fna", "ERS235603_16.fna", 
           "ERS396297_11.fna", "ERS396519_11.fna", "ERS473255_26.fna", 
           "ERS537218_9.fna", "ERS537235_19.fna", "ERS537328_30.fna", 
           "ERS537353_12.fna", "ERS608524_37.fna", "ERS608576_22.fna", 
           "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna", 
           "GCF_000508885.1_ASM50888v1_genomic.fna", 
           "GCF_001405615.1_13414_6_47_genomic.fna", 
           "GCF_900036035.1_RGNV35913_genomic.fna", 
           "LeChatelierE_2013__MH0074__bin.19.fa", "LiJ_2014__O2.UC28-1__bin.61.fa",
           "LiSS_2016__FAT_DON_8-22-0-0__bin.28.fa", "LoombaR_2017__SID1050_bax__bin.11.fa",
           "NielsenHB_2014__MH0094__bin.44.fa", "QinJ_2012__CON-091__bin.20.fa",
           "SRR4305229_bin.5.fa", "SRR5127401_bin.3.fa", "SRR5558047_bin.10.fa",
           "SRR6028281_bin.3.fa", "SRS075078_49.fna", "SRS103987_37.fna", 
           "SRS104400_110.fna", "SRS143598_15.fna", "SRS1719112_8.fna", 
           "SRS1719498_9.fna", "SRS1719577_6.fna", "SRS1735506_4.fna", 
           "SRS1735645_19.fna", "SRS294916_20.fna", "SRS476209_42.fna", 
           "VatanenT_2016__G80445__bin.9.fa", "VogtmannE_2016__MMRS43563715ST-27-0-0__bin.70.fa",
           "XieH_2016__YSZC12003_37172__bin.63.fa", "ZeeviD_2015__PNP_Main_232__bin.27.fa"]

rule all:
    input:
        # SOURMASH COMPARE OUTPUTS:
        "outputs/comp/all_filt_permanova_cosine.csv",
        "outputs/comp/all_filt_permanova_jaccard.csv",
        "outputs/comp/study_plt_all_filt_jaccard.pdf",
        "outputs/comp/diagnosis_plt_all_filt_jaccard.pdf",
        "outputs/comp/study_plt_all_filt_cosine.pdf",
        "outputs/comp/diagnosis_plt_all_filt_cosine.pdf",
        # VARIABLE SELECTION OUTPUTS:
        "outputs/filt_sig_hashes/count_total_hashes.txt",
        expand("outputs/vita_rf/{study}_vita_rf.RDS", study = STUDY),
        expand("outputs/vita_rf/{study}_vita_vars.txt", study = STUDY),
        expand("outputs/vita_rf/{study}_ibd_filt.csv", study = STUDY),
        # OPTIMAL RF OUTPUTS:
        expand('outputs/optimal_rf/{study}_optimal_rf.RDS', study = STUDY),
        # VARIABLE CHARACTERIZATION OUTPUTS:
        expand("outputs/gather/{study}_vita_vars_refseq.csv", study = STUDY),
        expand("outputs/gather/{study}_vita_vars_genbank.csv", study = STUDY),
        expand("outputs/gather/{study}_vita_vars_all.csv", study = STUDY),
        "outputs/gather_matches_hash_map/hash_to_genome_map_gather_all.csv",
        "aggregated_checkpoints/finished_collect_gather_vita_vars_all_sig_matches_lca_classify.txt",
        "aggregated_checkpoints/finished_collect_gather_vita_vars_all_sig_matches_lca_summarize.txt",
        # SPACEGRAPHCATS OUTPUTS:
        #expand("outputs/nbhd_read_sigs/{library}/{gather_genome}.cdbg_ids.reads.sig", library = LIBRARIES, gather_genome = GATHER_GENOMES),
        #expand("outputs/nbhd_reads_corncob/{gather_genome}_sig_ccs.tsv", gather_genome = GATHER_GENOMES),
        # PANGENOME SIGS
        "outputs/gather_sgc_pangenome/hash_to_genome_map_at_least_5_studies_pangenome.csv",
        "outputs/gather_sgc_pangenome/at_least_5_studies_vita_vars_all.csv",
        expand("outputs/gather_sgc_pangenome/{study}_vita_vars_all.csv", study = STUDY),
        expand("outputs/gather_sgc_pangenome/{study}_vita_vars_pangenome.csv", study = STUDY)


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
    conda: 'envs/env.yml'
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
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule fastp_trimmed_reads:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: "outputs/fastp_abundtrim/{library}.abundtrim.fastp.json"
    conda: "envs/fastp.yml"
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
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
    '''

rule compute_signatures:
    input: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output: "outputs/sigs/{library}.sig"
    conda: 'envs/env.yml'
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
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name greater_than_one_count_hashes --filename {input} {input}
    '''

rule filter_signatures_to_greater_than_1_hashes:
    input:
        filt_sig = "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig",
        sigs = "outputs/sigs/{library}.sig"
    output: "outputs/filt_sigs/{library}_filt.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_filtered_sigs:
    input: "outputs/filt_sigs/{library}_filt.sig"
    output: "outputs/filt_sigs_named/{library}_filt_named.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.library}_filt
    '''

rule describe_filtered_sigs:
    input: expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES)
    output: "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash signature describe --csv {output} {input}
    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/filt_sigs_named/{library}_filt_named.sig"
    output: "outputs/filt_sigs_named_csv/{library}_filt_named.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/filt_sigs_named_csv/{library}_filt_named.csv", library = LIBRARIES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'envs/r.yml'
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
    conda: 'envs/rf.yml'
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
    conda: 'envs/rf.yml'
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
    conda: 'envs/tuneranger.yml'
    script: "scripts/tune_rf.R"


############################################
## Predictive hash characterization - gather
############################################

rule convert_vita_vars_to_sig:
    input: "outputs/vita_rf/{study}_vita_vars.txt"
    output: "outputs/vita_rf/{study}_vita_vars.sig"
    conda: "envs/sourmash.yml"
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
    conda: 'envs/sourmash.yml'
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
    conda: 'envs/sourmash.yml'
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
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

rule merge_vita_vars_sig_all:
    input: expand("outputs/vita_rf/{study}_vita_vars.sig", study = STUDY)
    output: "outputs/vita_rf/vita_vars_merged.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig merge -o {output} {input}
    '''

rule combine_gather_vita_vars_all:
    output: "outputs/gather/vita_vars_all.csv"
    input: expand("outputs/gather/{study}_vita_vars_all.csv", study = STUDY)
    run:
        import pandas as pd
        
        li = []
        for filename in input:
            df = pd.read_csv(str(filename), index_col=None, header=0)
            df["study"] = str(filename)
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(str(output))


checkpoint collect_gather_vita_vars_all_sig_matches:
    input:
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json",
        csv="outputs/gather/vita_vars_all.csv"
    output: directory("outputs/gather_matches/")
    run:
        from sourmash import signature
        import pandas as pd

        # load gather results
        df = pd.read_csv(input.csv)

        # for each row, determine which database the result came from
        for index, row in df.iterrows():
            if row["filename"] == "inputs/gather_databases/almeida-mags-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.almeida-mags-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/genbank-d2-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.genbank-d2-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/nayfach-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.nayfach-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/pasolli-mags-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.pasolli-mags-k31/" + row["md5"]
            # open the signature, parse its name, and write the signature out to a new
            # folder
            sigfp = open(sigfp, 'rt')
            sig = signature.load_one_signature(sigfp)
            out_sig = str(sig.name())
            out_sig = out_sig.split('/')[-1]
            out_sig = out_sig.split(" ")[0]
            out_sig = "outputs/gather_matches/" + out_sig + ".sig"
            with open(str(out_sig), 'wt') as fp:
                signature.save_signatures([sig], fp)      


def aggregate_collect_gather_vita_vars_all_sig_matches(wildcards):
    checkpoint_output = checkpoints.collect_gather_vita_vars_all_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches/{genome}.sig", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.sig")).genome)
    return file_names
    

rule create_hash_genome_map_gather_vita_vars_all:
    input:
        #genomes = "outputs/gather_matches/{genome}.sig",
        genomes = aggregate_collect_gather_vita_vars_all_sig_matches,
        vita_vars = "outputs/vita_rf/vita_vars_merged.sig"
    output: "outputs/gather_matches_hash_map/hash_to_genome_map_gather_all.csv"
    run:
        from sourmash import signature
        import pandas as pd
        
        sigs = input.genomes
        # read in all genome signatures that had gather 
        # matches for the var imp hashes create a dictionary, 
        # where the key is the genome and the values are the minhashes
        genome_dict = {}
        for sig in sigs:
            sigfp = open(sig, 'rt')
            siglist = list(signature.load_signatures(sigfp))
            loaded_sig = siglist[0] 
            mins = loaded_sig.minhash.get_mins() # Get the minhashes 
            genome_dict[sig] = mins

        # read in vita variables
        sigfp = open(str(input.vita_vars), 'rt')
        vita_vars = sig = signature.load_one_signature(sigfp)
        vita_vars = vita_vars.minhash.get_mins() 

        # generate a list of all minhashes from all genomes
        all_mins = []
        for file in sigs:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[0]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes 
                all_mins += mins

        # define a function where if a hash is a value, 
        # return all key for which it is a value
        def get_all_keys_if_value(dictionary, hash_query):
            genomes = list()
            for genome, v in dictionary.items():
                if hash_query in v:
                    genomes.append(genome)
            return genomes

        # create a dictionary where each vita_vars hash is a key, 
        # and values are the genome signatures in which that hash
        # is contained
        vita_hash_dict = {}
        for hashy in vita_vars:
            keys = get_all_keys_if_value(genome_dict, hashy)
            vita_hash_dict[hashy] = keys

        # transform this dictionary into a dataframe and format the info nicely
        df = pd.DataFrame(list(vita_hash_dict.values()), index = vita_hash_dict.keys())
        df = df.reset_index()
        df = pd.melt(df, id_vars=['index'], var_name= "drop", value_name='genome')
        # remove tmp col drop
        df = df.drop('drop', 1)
        # drop duplicate rows in the df
        df = df.drop_duplicates()
        # write the dataframe to csv
        df.to_csv(str(output), index = False) 


rule download_sourmash_lca_db:
    output: "inputs/gather_databases/gtdb-release89-k31.lca.json.gz"
    shell:'''
    wget -O {output} https://osf.io/gs29b/download
    '''

rule sourmash_lca_classify_vita_vars_all_sig_matches:
    input:
        db = "inputs/gather_databases/gtdb-release89-k31.lca.json.gz",
        genomes = "outputs/gather_matches/{genome}.sig"
    output: "outputs/gather_matches_lca_classify/{genome}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash lca classify --db {input.db} --query {input.genomes} -o {output}
    '''

def aggregate_collect_gather_vita_vars_all_sig_matches_lca_classify(wildcards):
    checkpoint_output = checkpoints.collect_gather_vita_vars_all_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches_lca_classify/{genome}.csv", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.sig")).genome)
    return file_names

    
rule finished_collect_gather_vita_vars_all_sig_matches_lca_classify:
    input: aggregate_collect_gather_vita_vars_all_sig_matches_lca_classify
    output: "aggregated_checkpoints/finished_collect_gather_vita_vars_all_sig_matches_lca_classify.txt"
    shell:'''
    touch {output}
    '''

rule sourmash_lca_summarize_vita_vars_all_sig_matches:
    input:
        db = "inputs/gather_databases/gtdb-release89-k31.lca.json.gz",
        genomes = "outputs/gather_matches/{genome}.sig"
    output: "outputs/gather_matches_lca_summarize/{genome}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash lca summarize --db {input.db} --query {input.genomes} -o {output}
    '''

def aggregate_collect_gather_vita_vars_all_sig_matches_lca_summarize(wildcards):
    checkpoint_output = checkpoints.collect_gather_vita_vars_all_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches_lca_summarize/{genome}.csv", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.sig")).genome)
    return file_names

    
rule finished_collect_gather_vita_vars_all_sig_matches_lca_summarize:
    input: aggregate_collect_gather_vita_vars_all_sig_matches_lca_summarize
    output: "aggregated_checkpoints/finished_collect_gather_vita_vars_all_sig_matches_lca_summarize.txt"
    shell:'''
    touch {output}
    '''

###################################################
# Predictive hash characterization -- shared hashes
###################################################

rule at_least_5_of_6_hashes:
    """
    R script that takes as input the output vita vars,
    intersects the vars, and writes out to text file
    """
    input: expand("outputs/vita_rf/{study}_vita_vars.txt", study = STUDY)
    output: at_least_5 = "outputs/vita_rf/at_least_5_studies_vita_vars.txt"
    conda: 'envs/ggplot.yml'
    script: 'scripts/at_least_5_studies.R'


rule at_least_5_of_6_sig:
   """
   convert output of at_least_5_of_6_hashes to signature
   """
   input: "outputs/vita_rf/at_least_5_studies_vita_vars.txt"
   output: "outputs/vita_rf/at_least_5_studies_vita_vars.sig"
   conda: "envs/sourmash.yml"
   shell:'''
   python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name at_least_5_models --filename {input} {input}
   '''

rule at_least_5_of_6_gather:
    """
    run gather on the signature that contains hashes from
    at least 5 of 6 models
    """
    input:
        sig="outputs/vita_rf/at_least_5_studies_vita_vars.sig",
        db1="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db2="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db3="inputs/gather_databases/pasolli-mags-k31.sbt.json",
        db4="inputs/gather_databases/nayfach-k31.sbt.json",
    output: 
        csv="outputs/gather/at_least_5_studies_vita_vars.csv",
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --scaled 2000 -k 31 {input.sig} {input.db1} {input.db2} {input.db3} {input.db4}
    '''

checkpoint collect_gather_at_least_5_of_6_sig_matches:
    input:
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json",
        csv="outputs/gather/at_least_5_studies_vita_vars.csv"
    output: directory("outputs/gather_matches_loso_sigs/")
    run:
        from sourmash import signature
        import pandas as pd

        # load gather results
        df = pd.read_csv(input.csv)

        # for each row, determine which database the result came from
        for index, row in df.iterrows():
            if row["filename"] == "inputs/gather_databases/almeida-mags-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.almeida-mags-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/genbank-d2-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.genbank-d2-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/nayfach-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.nayfach-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/pasolli-mags-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.pasolli-mags-k31/" + row["md5"]
            # open the signature, parse its name, and write the signature out to a new
            # folder
            sigfp = open(sigfp, 'rt')
            sig = signature.load_one_signature(sigfp)
            out_sig = str(sig.name())
            out_sig = out_sig.split('/')[-1]
            out_sig = out_sig.split(" ")[0]
            out_sig = "outputs/gather_matches_loso_sigs/" + out_sig + ".sig"
            with open(str(out_sig), 'wt') as fp:
                signature.save_signatures([sig], fp)      


def aggregate_collect_gather_at_least_5_of_6_sig_matches(wildcards):
    checkpoint_output = checkpoints.collect_gather_at_least_5_of_6_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches_loso_sigs/{genome41}.sig", 
                        genome41 = glob_wildcards(os.path.join(checkpoint_output, "{genome41}.sig")).genome41)
    return file_names
    
rule compare_at_least_5_of_6_sigs:
    input: aggregate_collect_gather_at_least_5_of_6_sig_matches
    output: "outputs/comp_loso/comp_jaccard"
    conda: "envs/sourmash.yml"
    shell:''' 
    sourmash compare --ignore-abundance -k 31 -o {output} {input}
    '''

rule plot_at_least_5_of_6_sigs:
    input: "outputs/comp_loso/comp_jaccard"
    output: "outputs/comp_loso/comp_jaccard.matrix.pdf"
    params: out_dir = "outputs/comp_loso"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash plot --pdf --labels --output-dir {params.out_dir} {input} 
    '''

rule create_hash_genome_map_at_least_5_of_6_vita_vars:
    input:
        sigs = aggregate_collect_gather_at_least_5_of_6_sig_matches,
        gather = "outputs/gather/at_least_5_studies_vita_vars.csv",
    output: "outputs/gather_matches_loso_hash_map/hash_to_genome_map_at_least_5_studies.csv"
    run:
        files = input.sigs
         # load in all signatures that had gather matches and generate a list of all hashes 
        all_mins = []
        for file in files:
            sigfp = open(file, 'rt')
            siglist = list(signature.load_signatures(sigfp))
            loaded_sig = siglist[0] # sigs are from the sbt, only one sig in each (k=31)
            mins = loaded_sig.minhash.get_mins() # Get the minhashes 
            all_mins += mins # load all minhashes into a list

        # make all_mins a set
        all_mins = set(all_mins)

        # read in the gather matches as a dataframe
        gather_matches = pd.read_csv(input.gather)

        # loop through gather matches. For each match, read in it's
        # signature and generate an intersection of its hashes and all of 
        # the hashes that matched to gather. 
        # Then, remove those hashes from the all_mins list, and repeat
        sig_dict = {}
        all_sig_df = pd.DataFrame(columns=['level_0', '0'])
        for index, row in gather_matches.iterrows():
            sig = row["name"]
            # edit the name to match the signature names
            sig = str(sig)
            sig = sig.split('/')[-1]
            sig = sig.split(" ")[0]
            sig = sig + ".sig"
            sig = "outputs/gather_matches_loso_sigs/" + sig
            # load in signature
            sigfp = open(sig, 'rt')
            siglist = list(signature.load_signatures(sigfp))
            loaded_sig = siglist[0] # sigs are from the sbt, only one sig in each (k=31)
            mins = loaded_sig.minhash.get_mins() # Get the minhashes 
            mins = set(mins)
            # intersect all_mins list with mins from current signature
            intersect_mins = mins.intersection(all_mins)
            # add hashes owned by the signature to a dictionary 
            sig_dict[sig] = intersect_mins
            # convert into a dataframe
            sig_df = pd.DataFrame.from_dict(sig_dict,'index').stack().reset_index(level=0)
            # combine dfs
            all_sig_df = pd.concat([all_sig_df, sig_df], sort = False)
            # subtract intersect_mins from all_mins
            all_mins = all_mins - mins

        all_sig_df.columns = ['sig', 'tmp', 'hash']
        all_sig_df = all_sig_df.drop(['tmp'], axis = 1)
        all_sig_df = all_sig_df.drop_duplicates(keep = "first")
        all_sig_df.to_csv(str(output))


#############################################
# Spacegraphcats Genome Queries
#############################################

rule download_gather_match_genomes:
    output: "outputs/gather/gather_genomes_loso.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/tsu9c/download
    '''

#checkpoint untar_gather_match_genomes:
#    output:  directory("outputs/gather_matches_loso")
#    input:"outputs/gather/gather_genomes_loso.tar.gz"
#    params: outdir = "outputs/"
#    shell:'''
#    mkdir -p {params.outdir}
#    tar xf {input} -C {params.outdir}
#    '''

#def aggregate_untar_gather_match_genomes(wildcards):
#    # checkpoint_output produces the output dir from the checkpoint rule.
#    checkpoint_output = checkpoints.untar_gather_match_genomes.get(**wildcards).output[0]    
#    file_names = expand("outputs/gather_matches_loso/{gather_genome}.gz",
#                        gather_genome = glob_wildcards(os.path.join(checkpoint_output, "{gather_genome}.gz")).gather_genome)
#    return file_names

#rule spacegraphcats_gather_matches:
#    input: 
#        query = aggregate_untar_gather_match_genomes,
#        conf = "inputs/sgc_conf/{library}_r1_conf.yml",
#        reads = "outputs/abundtrim/{library}.abundtrim.fq.gz"
#    output:
#        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.fa.gz",
#        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.contigs.sig"
#    params: outdir = "outputs/sgc_genome_queries"
#    conda: "envs/spacegraphcats.yml"
#    shell:'''
#    spacegraphcats {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir}  
#    '''

rule untar_gather_match_genomes:
    output:  expand("outputs/gather_matches_loso/{gather_genome}.gz", gather_genome = GATHER_GENOMES)
    input:"outputs/gather/gather_genomes_loso.tar.gz"
    params: outdir = "outputs/"
    shell:'''
    mkdir -p {params.outdir}
    tar xf {input} -C {params.outdir}
    '''

rule spacegraphcats_gather_matches:
    input: 
        query = "outputs/gather_matches_loso/{gather_genome}.gz", 
        conf = "inputs/sgc_conf/{library}_r1_conf.yml",
        reads = "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output:
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.fa.gz",
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.contigs.sig"
    params: outdir = "outputs/sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    shell:'''
    spacegraphcats {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir}  
    '''
##############################################
## Pangenome signature/variable importance
##############################################

# use default contig sigs output by sgc to start. 

#rule calc_sig_nbhd_reads:
#    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.fa.gz"
#    output: "outputs/nbhd_read_sigs/{library}/{gather_genome}.cdbg_ids.reads.sig"
#    conda: "envs/sourmash.yml"
#    shell:'''
#    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} --merge {wildcards.library}_{wildcards.gather_genome} {input}
#    '''

rule merge_sgc_sigs_to_pangenome:
    input: expand("outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{{gather_genome}}.gz.contigs.sig", library = LIBRARIES)
    output: "outputs/sgc_pangenome_sigs/{gather_genome}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash signature merge -o {output} -k 31 {input}
    '''

rule rename_sgc_pangenome_sigs:
    input: "outputs/sgc_pangenome_sigs/{gather_genome}.sig"
    output: "outputs/sgc_pangenome_sigs/{gather_genome}_renamed.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig rename -o {output} {input} {wildcards.gather_genome}_pangenome
    '''

rule compare_sgc_pangenome_sigs:
    input: expand("outputs/sgc_pangenome_sigs/{gather_genome}_renamed.sig", gather_genome = GATHER_GENOMES)
    output: "outputs/sgc_pangenome_compare/pangenome_compare.comp"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compare -o {output} --ignore-abundance {input}
    '''

rule plot_sgc_pangenome_sigs:
    input: "outputs/sgc_pangenome_compare/pangenome_compare.comp"
    output: "outputs/sgc_pangenome_compare/pangenome_compare.comp.matrix.pdf"
    conda:"envs/sourmash.yml"
    shell:'''
    sourmash plot --labels {input}
    '''

rule compare_sgc_pangenome_sigs_csv:
    input: expand("outputs/sgc_pangenome_sigs/{gather_genome}_renamed.sig", gather_genome = GATHER_GENOMES)
    output: "outputs/sgc_pangenome_compare/pangenome_compare.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compare --csv {output} --ignore-abundance {input}
    '''

rule index_sgc_pangenome_sigs:
    input: expand("outputs/sgc_pangenome_sigs/{gather_genome}_renamed.sig", gather_genome = GATHER_GENOMES)
    output: "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash index -k 31 {output} {input}
    '''

rule gather_vita_vars_study_against_sgc_pangenome_sigs:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db = "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json"
    output: 
        csv="outputs/gather_sgc_pangenome/{study}_vita_vars_pangenome.csv",
        un="outputs/gather_sgc_pangenome/{study}_vita_vars_pangenome.un"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

rule gather_against_sgc_pangenome_sigs_plus_all_dbs:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db = "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json",
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json"
    output: 
        csv="outputs/gather_sgc_pangenome/{study}_vita_vars_all.csv",
        un="outputs/gather_sgc_pangenome/{study}_vita_vars_all.un"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db} {input.db1} {input.db2} {input.db3} {input.db4} 
    '''

rule at_least_5_of_6_gather_pangenome_sigs_plus_all_dbs:
    """
    run gather on the signature that contains hashes from
    at least 5 of 6 models
    """
    input:
        sig="outputs/vita_rf/at_least_5_studies_vita_vars.sig",
        db = "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json",
        db1="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db2="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db3="inputs/gather_databases/pasolli-mags-k31.sbt.json",
        db4="inputs/gather_databases/nayfach-k31.sbt.json",
    output: 
        csv="outputs/gather_sgc_pangenome/at_least_5_studies_vita_vars_all.csv",
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --scaled 2000 -k 31 {input.sig} {input.db} {input.db1} {input.db2} {input.db3} {input.db4}
    '''

rule at_least_5_of_6_gather_pangenome_sigs:
    """
    run gather on the signature that contains hashes from
    at least 5 of 6 models
    """
    input:
        sig="outputs/vita_rf/at_least_5_studies_vita_vars.sig",
        db = "outputs/sgc_pangenome_db/merged_sgc_sig.sbt.json"
    output: 
        csv="outputs/gather_sgc_pangenome/at_least_5_studies_vita_vars_pangenome.csv",
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

rule create_hash_genome_map_at_least_5_of_6_vita_vars_pangenome:
    input:
        sigs = expand("outputs/sgc_pangenome_sigs/{gather_genome}_renamed.sig", gather_genome = GATHER_GENOMES),
        gather="outputs/gather_sgc_pangenome/at_least_5_studies_vita_vars_pangenome.csv",
    output: "outputs/gather_sgc_pangenome/hash_to_genome_map_at_least_5_studies_pangenome.csv"
    run:
        files = input.sigs
         # load in all signatures that had gather matches and generate a list of all hashes 
        all_mins = []
        for file in files:
            sigfp = open(file, 'rt')
            siglist = list(signature.load_signatures(sigfp))
            loaded_sig = siglist[0] # sigs are from the sbt, only one sig in each (k=31)
            mins = loaded_sig.minhash.get_mins() # Get the minhashes 
            all_mins += mins # load all minhashes into a list

        # make all_mins a set
        all_mins = set(all_mins)

        # read in the gather matches as a dataframe
        gather_matches = pd.read_csv(input.gather)

        # loop through pangenome sigs in the order that they occurred as
        # gather matches. For each match, read in it's
        # signature and generate an intersection of its hashes and all of 
        # the hashes that matched to gather. 
        # Then, remove those hashes from the all_mins list, and repeat
        sig_dict = {}
        all_sig_df = pd.DataFrame(columns=['level_0', '0'])
        for index, row in gather_matches.iterrows():
            sig = row["name"]
            # edit the name to match the signature names
            sig = str(sig)
            sig = sig.split('/')[-1]
            sig = sig.split(" ")[0]
            sig = sig + ".sig"
            sig = "outputs/gather_sgc_pangenome/" + sig
            # load in signature
            sigfp = open(sig, 'rt')
            siglist = list(signature.load_signatures(sigfp))
            loaded_sig = siglist[0] # sigs are from the sbt, only one sig in each (k=31)
            mins = loaded_sig.minhash.get_mins() # Get the minhashes 
            mins = set(mins)
            # intersect all_mins list with mins from current signature
            intersect_mins = mins.intersection(all_mins)
            # add hashes owned by the signature to a dictionary 
            sig_dict[sig] = intersect_mins
            # convert into a dataframe
            sig_df = pd.DataFrame.from_dict(sig_dict,'index').stack().reset_index(level=0)
            # combine dfs
            all_sig_df = pd.concat([all_sig_df, sig_df], sort = False)
            # subtract intersect_mins from all_mins
            all_mins = all_mins - mins

        all_sig_df.columns = ['sig', 'tmp', 'hash']
        all_sig_df = all_sig_df.drop(['tmp'], axis = 1)
        all_sig_df = all_sig_df.drop_duplicates(keep = "first")
        all_sig_df.to_csv(str(output))


##############################################
## Pangenome differential abundance analysis
##############################################

rule diginorm_nbhd_reads:
    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_reads_diginorm/{library}/{gather_genome}.cdgb_ids.reads.diginorm.fa.gz"
    conda: "envs/env.yml"
    shell:'''
    normalize-by-median.py -k 20 -C 20 -M 16e9 --gzip -o {output} {input}
    '''

rule cat_diginorm_nbhd_reads:
    input: expand("outputs/nbhd_reads_diginorm/{library}/{{gather_genome}}.cdgb_ids.reads.diginorm.fa.gz", library = LIBRARIES)
    output: "outputs/nbhd_reads_diginorm_cat/{gather_genome}.cdbg_ids.reads.diginorm.fa.gz"
    shell:'''
    cat {input} > {output}
    '''

rule plass_nbhd_reads:
    input: "outputs/nbhd_reads_diginorm_cat/{gather_genome}.cdbg_ids.reads.diginorm.fa.gz"
    output: "outputs/nbhd_reads_plass/{gather_genome}.cdbg_ids.reads.plass.faa"
    conda: "envs/plass.yml"
    shell:'''
    plass assemble --threads 4 --min-length 25 {input} {output} tmp
    '''

rule plass_remove_stops:
    input: "outputs/nbhd_reads_plass/{gather_genome}.cdbg_ids.reads.plass.faa"
    output: "outputs/nbhd_reads_plass/{gather_genome}.cdbg_ids.reads.plass.faa.nostop.fa"
    conda: "envs/screed.yml"
    shell:'''
    python scripts/remove-stop-plass.py {input}
    '''

rule cdhit_plass:
    input: "outputs/nbhd_reads_plass/{gather_genome}.cdbg_ids.reads.plass.faa.nostop.fa"
    output: "outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "envs/plass.yml"
    shell:'''
    cd-hit -M 15500 -i {input} -o {output} -c .9
    '''

rule paladin_index_plass:
    input: "outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    output: "outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa.bwt"
    conda: "envs/plass.yml"
    shell: '''
    paladin index -r3 {input}
    '''

rule paladin_align_plass:
    input:
        indx="outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa.bwt",
        reads="outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.gz.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_reads_paladin/{library}/{gather_genome}.sam"
    params: indx = "outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "envs/plass.yml"
    shell:'''
    paladin align -f 125 -t 2 {params.indx} {input.reads} > {output}
    '''

rule samtools_view_paladin:
    output: "outputs/nbhd_reads_paladin/{library}/{gather_genome}.bam"
    input: "outputs/nbhd_reads_paladin/{library}/{gather_genome}.sam"
    conda: "envs/plass.yml"
    shell:'''
    samtools view -b {input} > {output}
    '''

rule salmon_paladin:
    output: "outputs/nbhd_reads_salmon/{library}/{gather_genome}_quant/quant.sf"
    input:
        cdhit="outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa",
        bam="outputs/nbhd_reads_paladin/{library}/{gather_genome}.bam"
    params:
        out="outputs/nbhd_reads_salmon/{library}/{gather_genome}_quant"
    conda: "envs/plass.yml"
    shell:'''
    salmon quant -t {input.cdhit} -l A -a {input.bam} -o {params.out} --minAssignedFrags 1 || touch {output}
    '''

rule install_corncob:
    output:
        corncob = "outputs/nbhd_reads_corncob/corncob_install.txt"
    conda: 'envs/corncob.yml'
    script: "scripts/install_corncob.R"

rule make_salmon_counts:
    input:
        quant= expand("outputs/nbhd_reads_salmon/{library}/{{gather_genome}}_quant/quant.sf", library = LIBRARIES),
        info = "inputs/working_metadata.tsv",
        mqc_fastp = "outputs/fastp_abundtrim/multiqc_data/mqc_fastp_filtered_reads_plot_1.txt",
    output:
        counts = "outputs/nbhd_reads_tximport/{gather_genome}_counts_raw.tsv",
    params: gather_genome = lambda wildcards: wildcards.gather_genome
    conda: 'envs/corncob.yml'
    script: "scripts/run_tximport.R"

rule filt_salmon_counts:
    input:
        counts = "outputs/nbhd_reads_tximport/{gather_genome}_counts_raw.tsv",
        info = "inputs/working_metadata.tsv",
        mqc_fastp = "outputs/fastp_abundtrim/multiqc_data/mqc_fastp_filtered_reads_plot_1.txt",
    output: 
        dim = "outputs/nbhd_reads_corncob/{gather_genome}_dim.tsv",
        filt = "outputs/nbhd_reads_tximport/{gather_genome}_counts.tsv"
    conda: 'envs/corncob.yml'
    script: "scripts/run_filt_salmon.R"

rule corncob_salmon:
    input:
        filt = "outputs/nbhd_reads_tximport/{gather_genome}_counts.tsv",
        info = "inputs/working_metadata.tsv",
        mqc_fastp = "outputs/fastp_abundtrim/multiqc_data/mqc_fastp_filtered_reads_plot_1.txt",
        corncob = "outputs/nbhd_reads_corncob/corncob_install.txt"
    output:
        all_ccs = "outputs/nbhd_reads_corncob/{gather_genome}_all_ccs.tsv",
        sig_ccs = "outputs/nbhd_reads_corncob/{gather_genome}_sig_ccs.tsv"
    conda: 'envs/corncob.yml'
    script: "scripts/run_corncob.R"

rule get_corncob_significant_aaseq_names:
    input: sig = "outputs/nbhd_reads_corncob/{gather_genome}_sig_ccs.tsv"
    output: names = "outputs/nbhd_reads_corncob/{gather_genome}_sig_ccs_names.txt"
    conda: 'envs/corncob.yml'
    script: "scripts/get_corncob_names.R"

rule grab_corncob_significant_aa_seqs:
    output: sig_faa = "outputs/nbhd_reads_corncob/{gather_genome}_sig_ccs.faa" 
    input:
        names = "outputs/nbhd_reads_corncob/{gather_genome}_sig_ccs_names.txt",
        fasta = "outputs/nbhd_reads_cdhit/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "envs/sourmash.yml"
    shell: '''
    scripts/extract-aaseq-matches.py {input.names} {input.fasta} > {output}
    '''

########################################
## PCoA
########################################

rule compare_signatures_cosine:
    input: 
        expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES),
    output: "outputs/comp/all_filt_comp_cosine.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compare -k 31 -p 8 --csv {output} {input}
    '''

rule compare_signatures_jaccard:
    input: 
        expand("outputs/filt_sigs_named/{library}_filt_named.sig", library = LIBRARIES),
    output: "outputs/comp/all_filt_comp_jaccard.csv"
    conda: "envs/sourmash.yml"
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
    script: "scripts/run_permanova.R"

rule permanova_cosine:
    input: 
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = "inputs/working_metadata.tsv",
        sig_info = "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    output: 
        perm = "outputs/comp/all_filt_permanova_cosine.csv"
    conda: "envs/vegan.yml"
    script: "scripts/run_permanova.R"

rule plot_comp_jaccard:
    input:
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        study = "outputs/comp/study_plt_all_filt_jaccard.pdf",
        diagnosis = "outputs/comp/diagnosis_plt_all_filt_jaccard.pdf"
    conda: "envs/ggplot.yml"
    script: "scripts/plot_comp.R"

rule plot_comp_cosine:
    input:
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = "inputs/working_metadata.tsv"
    output: 
        study = "outputs/comp/study_plt_all_filt_cosine.pdf",
        diagnosis = "outputs/comp/diagnosis_plt_all_filt_cosine.pdf"
    conda: "envs/ggplot.yml"
    script: "scripts/plot_comp.R"
