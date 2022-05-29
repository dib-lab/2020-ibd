# This snakefile was used to annotate the metapangenome species graphs.
# It requires the main Snakefile to be run in order to work, as it relies on outputs from that workflow (mainly "outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz", but it will also reuse the metapangenome species catlas if that already exists)
# This ended up being a simplified procedure to do annotation, as I directly specify the species and accessions I'm interested in exploring, and bind them together in a single wildcard.
# The class Checkpoint_GrabAccessions is still necessary to retrieve all of the genomes in GTDB for a given species.
# The other benefit is that this allowed me to specify that I was only interested in annotating 5 catlases.
# It is possible to automate this within the framework of the previous snakefile, but this was just simpler to run as a small separate piece.

import csv
import pandas as pd

TMPDIR = "/scratch/tereiter"

metadata = pd.read_csv("inputs/anno_metadata.tsv", sep = "\t", header = 0)
GTDB_SPECIES = metadata['species_no_space'].tolist()

outdir = "outputs"

ACC = ['GCF_008121495.1', 'GCF_900113155.1', 'GCF_000424325.1', 'GCF_002234575.2', 'GCF_005845215.1']
# make a variable that binds GTDB species to the accession that the queries were seeded with;
# catlases are named after the accession, so these need to be bound so that each catlas is only annotated by the correct species 
ACC_SPECIES = ['GCF_008121495.1--s__Ruminococcus_B-gnavus', 'GCF_900113155.1--s__Enterocloster-clostridioformis',
               'GCF_000424325.1--s__Enterocloster-clostridioformis_A', 'GCF_002234575.2--s__Enterocloster-bolteae',
               'GCF_005845215.1--s__Enterocloster-sp005845215']
           
class Checkpoint_GrabAccessions:
    """
    Define a class to simplify file specification from checkpoint 
    (e.g. solve for {acc} wildcard without needing to specify a function
    for each arm of the DAG that uses the acc wildcard). 
    This approach is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self, gtdb_species):
        acc_csv = f'outputs/genbank/{gtdb_species}.x.genbank.gather.csv'
        assert os.path.exists(acc_csv)

        genome_accs = []
        with open(acc_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['name']
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {acc_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of rule 'grab_species_accessions'; 
        # this will trigger exception until that rule has been run.
        checkpoints.grab_species_accessions.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(w.gtdb_species)

        p = expand(self.pattern, acc=genome_accs, **w)
        return p


rule all:
    input:
        expand("outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc_species}_distinct_enriched.tsv" , acc_species = ACC_SPECIES),
        expand("outputs/sgc_pangenome_catlases_selected_marker_gene_bandage/{acc}_done_marker_genes.txt", acc = ACC),
        expand("outputs/sgc_pangenome_catlases_selected_marker_gene_bwa/{acc}_done_marker_genes.txt", acc = ACC)

#################################################################
## Generate pangenome using reference sequences
#################################################################

checkpoint grab_species_accessions:
    input:
        lineages="/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.csv",
        metadata="inputs/anno_metadata.tsv"
    output: csv="outputs/genbank/{gtdb_species}.x.genbank.gather.csv",
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    script: "scripts/grab_species_accessions.R"

rule make_genome_info_csv:
    output:
        csvfile = 'outputs/genbank_genomes/{acc}.info.csv'
    conda: "envs/genbank_genomes.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell: """
        python scripts/genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
    """

rule download_matching_genome_wc:
    input:
        csvfile = ancient('outputs/genbank_genomes/{acc}.info.csv')
    output:
        genome = "outputs/genbank_genomes/{acc}_genomic.fna.gz"
    resources:
        mem_mb = 500
    threads: 1
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading genome for acc {acc}/{name} from NCBI...", file=sys.stderr)
            with open(output.genome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.genome}",
                        file=sys.stderr)

rule bakta_download_db:
    #output: "inputs/bakta_db/db/version.json"
    output: "/home/tereiter/github/2022-microberna/inputs/bakta_db/db/version.json" # note download already exists on cluster; change uncomment previous line and delete this line if download doesn't already exist. 
    threads: 1
    resources: 
        mem_mb = 4000,
        time_min = 240
    params: outdir = "inputs/bakta_db"
    conda: "envs/bakta.yml"
    shell:'''
    bakta_db download --output {params.outdir}
    '''
                       
rule bakta_annotate_gtdb_genomes:
    input: 
        fna=ancient("outputs/genbank_genomes/{acc}_genomic.fna.gz"),
        db="/home/tereiter/github/2022-microberna/inputs/bakta_db/db/version.json",
    output: 
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.faa",
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.gff3",
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna",
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.ffn",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        time_min = 40
    benchmark: "benchmarks/gtdb_genomes/bakta_{gtdb_species}/{acc}.txt"
    conda: 'envs/bakta.yml'
    params: 
        dbdir="/home/tereiter/github/2022-microberna/inputs/bakta_db/db/",
        outdirb = lambda wildcards: f'{outdir}/gtdb_genomes_bakta/' + wildcards.gtdb_species,
    threads: 1
    shell:'''
    bakta --threads {threads} --db {params.dbdir} --prefix {wildcards.acc} --output {params.outdirb} \
        --locus {wildcards.acc} --locus-tag {wildcards.acc} --keep-contig-headers {input.fna}
    '''

rule cat_annotated_sequences:
    input: ancient(Checkpoint_GrabAccessions(outdir + "/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.ffn"))
    output: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_all_annotated_seqs.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min = 20
    benchmark:"benchmarks/gtdb_genomes/cat_annotated_seqs_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''
    
rule cluster_annotated_sequences:
    input: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_all_annotated_seqs.fa"
    output: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa"
    threads: 1
    resources: 
        mem_mb=16000,
        time_min=30
    conda: "envs/cdhit.yml"
    benchmark:"benchmarks/gtdb_genomes/cluster_annotated_{gtdb_species}.txt"
    shell:'''
    cd-hit-est -c .95 -d 0 -i {input} -o {output}
    '''

####################################################################
## optional: ortholog annotation (else only bakta annots used)
####################################################################

rule translate_clustered_sequences_for_annotations:
    input: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa" 
    output: outdir + '/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.faa'
    conda: 'envs/emboss.yml'
    resources:
        mem_mb = 16000,
        time_min=20
    benchmark: "benchmarks/gtdb_genomes/transeq_{gtdb_species}.txt"
    threads: 2
    shell:'''
    transeq {input} {output}
    '''

rule eggnog_download_db:
    output: "/home/tereiter/github/2022-microberna/inputs/eggnog_db/eggnog.db"
    threads: 1   
    resources: 
        mem_mb = 4000,
        time_min=2880
    params: datadir = "inputs/eggnog_db"
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir {params.datadir}
    '''

rule eggnog_annotate_clustered_sequences:
    input: 
        faa = outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.faa",
        db = '/home/tereiter/github/2022-microberna/inputs/eggnog_db/eggnog.db'
    output: outdir + "/gtdb_genomes_annotated_comb_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations"
    conda: 'envs/eggnog.yml'
    resources:
        mem_mb = 32000,
        time_min=2880
    benchmark: "benchmarks/gtdb_genomes/eggnog_{gtdb_species}.txt"
    threads: 8
    params: 
        outdire = lambda wildcards: outdir + "/gtdb_genomes_annotated_comb_eggnog/" + wildcards.gtdb_species,
        dbdir = "/home/tereiter/github/2022-microberna/inputs/eggnog_db"
    shell:'''
    mkdir -p tmp/
    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.gtdb_species} \
       --output_dir {params.outdire} -m hmmer -d none --tax_scope auto \
       --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       -d 2 --data_dir {params.dbdir}
    '''

#####################################################################
## spacegraphcats multifastaqueries
#####################################################################

rule sketch_pangenome_reference:
    input: "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa",
    output: "outputs/gtdb_genomes_annotated_comb_sigs/{gtdb_species}_clustered_annotated_seqs.sig"
    resources:
        mem_mb = 500
    threads: 1
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=31,scaled=1000,abund -o {output} {input}
    '''
    
rule make_sgc_pangenome_multifasta_conf_files:
    input:
        reads = "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_genome_queries_hardtrim/{acc}.hardtrim.fa.gz",
        ref_genes = "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa",
        ref_sig = "outputs/gtdb_genomes_annotated_comb_sigs/{gtdb_species}_clustered_annotated_seqs.sig"
    output:
        conf = "outputs/sgc_conf/{acc}--{gtdb_species}_r10_multifasta_conf.yml"
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
multifasta_reference:
- {input.ref_genes}
multifasta_scaled: 1000
multifasta_query_sig: {input.ref_sig}
""", file=fp)


rule spacegraphcats_pangenome_catlas_multifasta_annotate:
    input:
        conf = "outputs/sgc_conf/{acc}--{gtdb_species}_r10_multifasta_conf.yml",
        catlas = "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
    output: 
        #annot="/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_annot.csv",
        #record="/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_by_record.csv",
        touch = "outputs/sgc_pangenome_catlases/{acc}--{gtdb_species}_done.txt"
    params: 
        outdir = "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/",
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 32000
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} multifasta_query --nolock --outdir {params.outdir} --rerun-incomplete
    touch {output.touch}
    '''

rule join_annotations_to_dominating_set_differential_abundance_analysis_results:
    input:
        cdbg_to_pieces="/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10/cdbg_to_pieces.csv",
        eggnog= "outputs/gtdb_genomes_annotated_comb_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations",
        cdbg_annot= "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_multifasta/multifasta.cdbg_annot.csv",
        sig_ccs = "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases_corncob/{acc}_sig_ccs.tsv",
        dom_info="/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_abund/dom_info.tsv"
    output:
        sig_ccs_annot="outputs/sgc_pangenome_catlases_corncob_annotations/{acc}--{gtdb_species}_sig_ccs_annotated_no_cdbg.tsv"
    threads: 1
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb=8000
    script: "scripts/join_annotations_to_dominating_set_differential_abundance_analysis_results.R"

rule summarize_function_and_enrichment_on_dominating_set_differential_abundance_analysis_results:
    input:
        sig_ccs_annot="outputs/sgc_pangenome_catlases_corncob_annotations/{acc}--{gtdb_species}_sig_ccs_annotated_no_cdbg.tsv",
        marker_kos = "inputs/kegg_scg_13059_2015_610_MOESM1_ESM_filtered.csv",
        ox_kos = "inputs/ox_stress.tsv",
    output: 
        enriched_distinct="outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc}--{gtdb_species}_distinct_enriched.tsv" , 
        enriched_overlap="outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc}--{gtdb_species}_overlapping_enriched.tsv",
        overlap_marker_kos="outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc}--{gtdb_species}_overlapping_marker_kos.tsv",
        overlap_df = "outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc}--{gtdb_species}_overlap_kos_df.tsv",
        distinct_ox = "outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc}--{gtdb_species}_distinct_ox_kos.tsv",
        distinct_df = "outputs/sgc_pangenome_catlases_corncob_annotation_analysis/{acc}--{gtdb_species}_distinct_kos_df.tsv"
    threads: 1
    resources:
        mem_mb=8000,
        time_min = 120
    conda: "envs/clusterprofiler.yml"
    script: "scripts/summarize_function_and_enrichment_on_dominating_set_differential_abundance_analysis_results.R"

checkpoint sgc_pangenome_catlases_extract_marker_gene_reads:
    input:
        conf = "outputs/sgc_conf/{acc}_r10_extract_marker_genes_conf.yml",
        catlas = "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10/catlas.csv"
    output: directory("/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_search_oh0")
    params: 
        outdir = "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/",
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 16000
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_reads --nolock --outdir {params.outdir} --rerun-incomplete -n
    '''

rule index_query:
    input: "outputs/sgc_pangenome_catlases_selected_marker_genes/{marker_gene}.fa"
    output: "outputs/sgc_pangenome_catlases_selected_marker_genes/{marker_gene}.fa.bwt"
    conda: "envs/bwa.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    bwa index {input}
    ''' 

rule map_reads_back_to_query:
    input:
        query="outputs/sgc_pangenome_catlases_selected_marker_genes/{marker_gene}.fa",
        query_bwt="outputs/sgc_pangenome_catlases_selected_marker_genes/{marker_gene}.fa.bwt",
        reads="/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_search_oh0/{marker_gene}.fa.cdbg_ids.reads.gz"
    output:"outputs/sgc_pangenome_catlases_selected_marker_gene_bwa/{acc}_r10/{marker_gene}.bam"
    conda: "envs/bwa.yml"
    threads: 4
    resources: mem_mb = 4000
    shell:'''
    bwa mem -p -t {threads} {input.query} {input.reads} | samtools sort -o {output} -
    '''

rule index_bam:
    input:"outputs/sgc_pangenome_catlases_selected_marker_gene_bwa/{acc}_r10/{marker_gene}.bam"
    output:"outputs/sgc_pangenome_catlases_selected_marker_gene_bwa/{acc}_r10/{marker_gene}.bam.bai"
    conda: "envs/bwa.yml"
    threads: 1
    resources: mem_mb = 4000
    shell:'''
    samtools index {input}
    '''

def checkpoint_sgc_pangenome_catlases_extract_marker_gene_reads_2(wildcards):
    # Expand checkpoint to get fasta names, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.sgc_pangenome_catlases_extract_marker_gene_reads.get(**wildcards).output[0]    
    file_names = expand("outputs/sgc_pangenome_catlases_selected_marker_gene_bwa/{{acc}}_r10/{marker_gene}.bam.bai",
                        marker_gene = glob_wildcards(os.path.join(checkpoint_output, "{marker_gene}.fa.cdbg_ids.reads.gz")).marker_gene)
    return file_names

rule dummy_solve_marker_gene_bwa:
    input: checkpoint_sgc_pangenome_catlases_extract_marker_gene_reads_2
    output: touch("outputs/sgc_pangenome_catlases_selected_marker_gene_bwa/{acc}_done_marker_genes.txt")


# Ran the following code and have results, but not the most useful visualization; use mapping back to query instead.
rule bcalm_marker_gene_nbhd:
    input: "/group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_pangenome_catlases/{acc}_k31_r10_search_oh0/{marker_gene}.fa.cdbg_ids.reads.gz"
    output: "outputs/sgc_pangenome_catlases_selected_marker_gene_bcalm/{acc}_r10/{marker_gene}.fa.cdbg_ids.reads.gz.unitigs.fa"
    params: prefix = "outputs/sgc_pangenome_catlases_selected_marker_gene_bcalm/{acc}_r10/{marker_gene}.fa.cdbg_ids.reads.gz"
    resources:
        mem_mb = 4000,
        time_min = 120
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    bcalm -in {input} -out-dir outputs/sgc_pangenome_catlases_selected_marker_gene_bcalm/{wildcards.acc}_r10 -kmer-size 31 -abundance-min 1 -out {params.prefix}
    '''

rule convert_to_gfa:
    input: "outputs/sgc_pangenome_catlases_selected_marker_gene_bcalm/{acc}_r10/{marker_gene}.fa.cdbg_ids.reads.gz.unitigs.fa"
    output: "outputs/sgc_pangenome_catlases_selected_marker_gene_bcalm/{acc}_r10/{marker_gene}.fa.cdbg_ids.reads.gz.unitigs.gfa"
    resources:
        mem_mb = 2000,
        time_min = 60
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    python ./scripts/convertToGFA.py {input} {output} 31
    '''

rule bandage_plot_gfa_marker_genes:
    input: 
        gfa="outputs/sgc_pangenome_catlases_selected_marker_gene_bcalm/{acc}_r10/{marker_gene}.fa.cdbg_ids.reads.gz.unitigs.gfa",
        blastquery="outputs/sgc_pangenome_catlases_selected_marker_genes/{marker_gene}.fa"
    output: "outputs/sgc_pangenome_catlases_selected_marker_gene_bandage/{acc}_r10/{marker_gene}.fa.cdbg_ids.reads.gz.unitigs.png"
    resources:
        mem_mb = 4000,
        time_min = 120
    threads: 1
    conda: "envs/bandage.yml"
    shell:'''
    XDG_RUNTIME_DIR=~/tmp # REQUIRES THAT THIS FOLDER EXISTS
    Bandage image {input.gfa} {output} --query {input.blastquery}    
    '''

def checkpoint_sgc_pangenome_catlases_extract_marker_gene_reads_1(wildcards):
    # Expand checkpoint to get fasta names, which will be used as queries for spacegraphcats extract
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.sgc_pangenome_catlases_extract_marker_gene_reads.get(**wildcards).output[0]    
    file_names = expand("outputs/sgc_pangenome_catlases_selected_marker_gene_bandage/{{acc}}_r10/{marker_gene}.fa.cdbg_ids.reads.gz.unitigs.png",
                        marker_gene = glob_wildcards(os.path.join(checkpoint_output, "{marker_gene}.fa.cdbg_ids.reads.gz")).marker_gene)
    return file_names

rule dummy_solve_marker_gene_blast:
    input: checkpoint_sgc_pangenome_catlases_extract_marker_gene_reads_1
    output: touch("outputs/sgc_pangenome_catlases_selected_marker_gene_bandage/{acc}_done_marker_genes.txt")


