# These rules weren't run for annotation of the metapangenome species graphs, but represent a work-in-progress set of rules that could be used to achieve the highest level of annotation on the catlas.
# I ended up using a more straightforward approach that only relied on genes from genomes in GTDB to cut down on the number of rules run.
# The strategy I used is encoded in annotate_metapangenome_species_graphs.snakefile.
# A precursor to the set of rules in this file is in roary.snakefile



# this class creates a mapping between acc:roary_acc,
# so it might just autosolve the connection between those wildcards
# and run the right things. If not, I guess try and make a dictionary
# in this class and use it somehow somewhere?
# dictionary would go in `with open(prefetch_csv` and init a dict 
# that would be available outside the class
class Checkpoint_RoaryPrefetchResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {roary_acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_roary_accs(self, acc): 
        prefetch_csv = f'outputs/roary_prefetch/{acc}_prefetch_filtered.csv'
        assert os.path.exists(prefetch_csv)

        roary_accs = []
        with open(prefetch_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               roary_acc = row['match_name'].split(' ')[0]
               roary_accs.append(roary_acc)
        print(f'loaded {len(roary_accs)} accessions from {prefetch_csv}.')

        return roary_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'touch_roary_filter_prefetch...';
        # this will trigger exception until that rule has been run.
        checkpoints.touch_roary_filter_prefetch_shared_assemblies_vs_refseq.get(**w)

        # parse accessions in gather output file
        roary_accs = self.get_roary_accs(**w) # w.accs has a hard time getting passed in here

        p = expand(self.pattern, roary_acc=roary_accs, **w)
        return p

####################################################
## Prepare multifasta reference annotation gene sets
####################################################

#rule sgc_genome_queries_repair:
#    input: "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{acc}_genomic.fna.gz.clean.fa.gz.cdbg_ids.reads.gz"
#    output:
#        pair="outputs/sgc_genome_queries_repaired/{library}_{acc}.repaired.fq.gz",
#        single="outputs/sgc_genome_queries_repaired/{library}_{acc}.orphaned.fq.gz"
#    threads: 1
#    resources: mem_mb = 4000
#    conda: "envs/bbmap.yml"
#    shell:'''
#    repair.sh in1={input} out={output.pair} outs={output.single} repair
#    '''
#    
#rule sgc_genome_queries_megahit:
#    input: "outputs/sgc_genome_queries_repaired/{library}_{acc}.repaired.fq.gz"
#    output: "outputs/sgc_genome_queries_megahit/{library}_{acc}.contigs.fa"
#    threads: 2
#    resources: mem_mb = 8000
#    conda: 'envs/megahit.yml'
#    shell:'''
#    megahit --12 {input} -t {threads} \
#        --out-dir outputs/tmp_megahit/{wildcards.library}_{wildcards.acc} \
#        --out-prefix {wildcards.library}_{wildcards.acc}
#    mv  outputs/tmp_megahit/{wildcards.library}_{wildcards.acc}/{wildcards.library}_{wildcards.acc}.contigs.fa {output}
#    rm -rf outputs/tmp_megahit/{wildcards.library}_{wildcards.acc}
#    '''
#    
#rule sgc_genome_queries_megahit_prokka:
#    output: 
#        ffn = 'outputs/sgc_genome_queries_megahit_prokka/{library}_{acc}.ffn',
#        faa = 'outputs/sgc_genome_queries_megahit_prokka/{library}_{acc}.faa',
#        gff = 'outputs/sgc_genome_queries_megahit_prokka/{library}_{acc}.gff',
#    input: 'outputs/sgc_genome_queries_megahit/{library}_{acc}.contigs.fa'
#    conda: 'envs/prokka.yml'
#    threads: 1
#    resources: mem_mb = 4000
#    params: 
#        output_folder = "outputs/sgc_genome_queries_megahit_prokka"
#    shell:'''
#    prokka {input} --outdir {params.output_folder} \
#       --prefix {wildcards.library}_{wildcards.acc} --metagenome --force \
#       --locustag {wildcards.library}_{wildcards.acc} --cpus {threads} # || touch {output.ffn}
    #touch {output.faa}
#    '''
#
## DOWNLOAD AND ANNOTATE ISOLATES FROM SPECIES IN SHARED ASSEMBLIES
#
#rule compute_signatures_shared_assemblies:
#    input: ancient('outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz')
#    output: 'outputs/charcoal/{acc}_genomic.fna.gz.clean.sig'
#    conda: "envs/sourmash.yml"
#    threads: 1
#    resources: mem_mb = 2000
#    shell:'''
#    sourmash compute -k 31 --scaled 2000 -o {output} {input}
#    '''
#
#rule roary_prefetch_shared_assemblies_vs_refseq:
#    input: 
#        query='outputs/charcoal/{acc}_genomic.fna.gz.clean.sig',
#        db = '/home/irber/sourmash_databases/outputs/sbt/refseq-bacteria-x1e6-k31.sbt.zip'
#    output: "outputs/roary_prefetch/{acc}_prefetch.csv"
#    conda: "envs/sourmash.yml"
#    threads: 1
#    resources: mem_mb = 16000
#    shell:'''
#    sourmash prefetch -k 31 -o {output} {input.query} {input.db} 
#    '''
#
#rule roary_filter_prefetch_shared_assemblies_vs_refseq:
## filter to jaccard of >= 0.1
#    input: prefetch="outputs/roary_prefetch/{acc}_prefetch.csv"
#    output: filt="outputs/roary_prefetch/{acc}_prefetch_filtered.csv"
#    conda: "envs/tidy.yml"
#    threads: 1
#    resources: mem_mb = 4000
#    script: "scripts/roary_filter_prefetch.R"
#
## input needs to be output from previous rule
## check that they all exist, and then touch a dummy file
#checkpoint touch_roary_filter_prefetch_shared_assemblies_vs_refseq:
#    input: "outputs/roary_prefetch/{acc}_prefetch_filtered.csv"
#    output: touch("outputs/roary_prefetch/.{acc}_roary_dummy.txt")
#    threads: 1
#    resources: mem_mb = 4000
#
#rule roary_combine_filtered_prefetch_results:
#    input: filt = Checkpoint_GatherResults("outputs/roary_prefetch/{acc}_prefetch_filtered.csv")
#    output: combined = "outputs/genbank/roary_prefetch.x.genbank.gather.csv"
#    conda: "envs/tidy.yml"
#    threads: 1
#    resources: mem_mb = 4000
#    script: "scripts/roary_combine.R"
#
#rule roary_make_conf_file_for_genome_grist:
#    input: "outputs/genbank/roary_prefetch.x.genbank.gather.csv"
#    output: conf="outputs/roary_genome_grist_conf/roary_genome_grist_conf.yml"
#    threads: 1
#    resources: mem_mb = 1000
#    run:
#        with open(output.conf, 'wt') as fp:
#           print(f"""\
#sample:
#- roary_prefetch
#outdir: outputs
#""", file=fp)
#
## in the Checkpoint_Roary class, we can create a dictionary of acc:roary_acc to
## handle stupid genome-grist dumbness. Then can use lambda wildcards to get 
## roary_acc values from acc key in dictionary. 
#rule roary_genome_grist:
#    input:
#        conf="outputs/roary_genome_grist_conf/roary_genome_grist_conf.yml"
#    output: 
#        genome = "outputs/genbank_genomes_roary/{roary_acc}_genomic.fna.gz"
#    conda: "envs/genome-grist.yml"
#    resources: mem_mb = 8000
#    threads: 1
#    shell:'''
#    genome-grist run {input.conf} --until make_sgc_conf --nolock
#    mv genbank_genomes/ outputs/genbank_genomes_roary
#    '''
#
#rule roary_prokka:
#    output: 
#        ffn = 'outputs/roary_prokka/{acc}/{roary_acc}.ffn',
#        faa = 'outputs/roary_prokka/{acc}/{roary_acc}.faa',
#        gff = "outputs/roary_prokka/{acc}/{roary_acc}.gff"
#    input: "outputs/genbank_genomes_roary/{roary_acc}_genomic.fna.gz"
#    conda: 'envs/prokka.yml'
#    resources: mem_mb = 8000
#    threads: 2
#    params: 
#        outdir = lambda wildcards: 'outputs/roary_prokka/' + wildcards.acc + "/",
#        prefix = lambda wildcards: wildcards.roary_acc,
#        gzip = lambda wildcards: "outputs/genbank_genomes_roary/" + wildcards.roary_acc + "_genomic.fna.gz"
#    shell:'''
#    gunzip {input}
#    prokka {params.gzip} --outdir {params.outdir} --prefix {params.prefix} --metagenome --force --locustag {params.prefix} --cpus {threads} --centre X --compliant
#    mv {params.prefix}.ffn {output.ffn}
#    mv {params.prefix}.faa {output.faa}
#    mv {params.prefix}.gff {output.gff}
#    gzip {params.gzip}
#    '''
#
#rule roary:
#    input: 
#        gff1 = Checkpoint_RoaryPrefetchResults('outputs/roary_prokka/{{acc}}/{roary_acc}.gff'),
#        dummy = "outputs/roary_prefetch/.{acc}_roary_dummy.txt",
#        gff2 = expand('outputs/sgc_genome_queries_megahit_prokka/{library}_{{acc}}.gff', library = LIBRARIES)
#    output: 'outputs/roary/{acc}/pan_genome_reference.fa'
#    conda: 'envs/roary.yml'
#    resources: mem_mb = 32000
#    threads: 4
#    params: outdir = lambda wildcards: "outputs/roary/" + wildcards.acc + "/"
#    shell:'''
#    roary -e -n -f {params.outdir} -p {threads} -z {input.gff1} {input.gff2}
#    '''
#
#rule roary_signature:
#    input: "outputs/roary/{acc}/pan_genome_reference.fa"
#    output: "outputs/roary/{acc}/pan_genome_reference.sig"
#    conda: 'envs/sourmash.yml'
#    resources: mem_mb = 4000
#    threads: 1
#    shell:'''
#    sourmash compute -k 31 --scaled 2000 -o {output} {input}
#    '''

#rule roary_translate:
#    input: 'outputs/roary/{acc}/pan_genome_reference.fa' 
#    output: 'outputs/roary/{acc}/pan_genome_reference.faa' 
#    conda: 'envs/emboss.yml'
#    resources: mem_mb = 16000
#    threads: 1
#    shell:'''
#    transeq {input} {output}
#    '''
#
#rule roary_eggnog:
#    input: 
#        faa = 'outputs/roary/{acc}/pan_genome_reference.faa',
#        db = 'inputs/eggnog_db/eggnog.db'
#    output: "outputs/roary_eggnog/{acc}.emapper.annotations"
#    conda: 'envs/eggnog.yml'
#    resources:
#        mem_mb = 64000
#    threads: 8
#    params: 
#        outdir = "outputs/roary_eggnog/",
#        dbdir = "inputs/eggnog_db"
#    shell:'''
#    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.acc} --output_dir {params.outdir} -m diamond -d none --tax_scope auto --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0 --override --temp_dir tmp/ -d bact --data_dir {params.dbdir}
#    '''
#
####################################################
## Query by multifasta eggnog gene annotation
####################################################

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
