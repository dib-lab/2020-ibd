# ================
import csv
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
#=========================
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

rule all:
    input: 
        Checkpoint_GatherResults("outputs/roary/{acc}/pan_genome_reference.fa")
#=======================================
checkpoint gather_gtdb_rep_to_shared_assemblies:
    output: 
        gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        gather_all_shared = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies_all.gather.csv"
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/gather_gtdb_rep_to_shared_assemblies.R"


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
    '''

#=========================================

rule compute_signatures_shared_assemblies:
    input: 'outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz'
    output: 'outputs/charcoal/{acc}_genomic.fna.gz.clean.sig'
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compute -k 31 --scaled 2000 -o {output} {input}
    '''

rule roary_prefetch_shared_assemblies_vs_refseq:
    input: 
        query='outputs/charcoal/{acc}_genomic.fna.gz.clean.sig',
        db = '/home/irber/sourmash_databases/outputs/sbt/refseq-bacteria-x1e6-k31.sbt.zip'
    output: "outputs/roary_prefetch/{acc}_prefetch.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash prefetch -k 31 -o {output} {input.query} {input.db} 
    '''

rule roary_filter_prefetch_shared_assemblies_vs_refseq:
# filter to jaccard of >= 0.1
    input: prefetch="outputs/roary_prefetch/{acc}_prefetch.csv"
    output: filt="outputs/roary_prefetch/{acc}_prefetch_filtered.csv"
    conda: "envs/tidy.yml"
    script: "scripts/roary_filter_prefetch.R"

# input needs to be output from previous rule
# check that they all exist, and then touch a dummy file
checkpoint touch_roary_filter_prefetch_shared_assemblies_vs_refseq:
    input: "outputs/roary_prefetch/{acc}_prefetch_filtered.csv"
    output: touch("outputs/roary_prefetch/.{acc}_roary_dummy.txt")

rule roary_combine_filtered_prefetch_results:
    input: filt = Checkpoint_GatherResults("outputs/roary_prefetch/{acc}_prefetch_filtered.csv")
    output: combined = "outputs/genbank/roary_prefetch.x.genbank.gather.csv"
    conda: "envs/tidy.yml"
    script: "scripts/roary_combine.R"

rule roary_make_conf_file_for_genome_grist:
    input: "outputs/genbank/roary_prefetch.x.genbank.gather.csv"
    output: conf="outputs/roary_genome_grist_conf/roary_genome_grist_conf.yml"
    run:
        with open(output.conf, 'wt') as fp:
           print(f"""\
sample:
- roary_prefetch
outdir: outputs
""", file=fp)

# in the Checkpoint_Roary class, we can create a dictionary of acc:roary_acc to
# handle stupid genome-grist dumbness. Then can use lambda wildcards to get 
# roary_acc values from acc key in dictionary. 
rule roary_genome_grist:
    input:
        conf="outputs/roary_genome_grist_conf/roary_genome_grist_conf.yml"
    output: 
        genome = "outputs/genbank_genomes_roary/{roary_acc}_genomic.fna.gz"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
    mv genbank_genomes/ outputs/genbank_genomes_roary
    '''

rule roary_prokka:
    output: 
        ffn = 'outputs/roary_prokka/{acc}/{roary_acc}.ffn',
        faa = 'outputs/roary_prokka/{acc}/{roary_acc}.faa',
        gff = "outputs/roary_prokka/{acc}/{roary_acc}.gff"
    input: "outputs/genbank_genomes_roary/{roary_acc}_genomic.fna.gz"
    conda: 'envs/prokka.yml'
    resources:
        mem_mb = 8000
    threads: 2
    params: 
        outdir = lambda wildcards: 'outputs/roary_prokka/' + wildcards.acc + "/",
        prefix = lambda wildcards: wildcards.roary_acc,
        gzip = lambda wildcards: "outputs/genbank_genomes_roary/" + wildcards.roary_acc + "_genomic.fna.gz"
    shell:'''
    gunzip {input}
    prokka {params.gzip} --outdir {params.outdir} --prefix {params.prefix} --metagenome --force --locustag {params.prefix} --cpus {threads} --centre X --compliant
    mv {params.prefix}.ffn {output.ffn}
    mv {params.prefix}.faa {output.faa}
    mv {params.prefix}.gff {output.gff}
    gzip {params.gzip}
    '''

rule roary:
    input: 
        gff = Checkpoint_RoaryPrefetchResults('outputs/roary_prokka/{acc}/{roary_acc}.gff'),
        dummy = "outputs/roary_prefetch/.{acc}_roary_dummy.txt"
    output: 'outputs/roary/{acc}/pan_genome_reference.fa'
    conda: 'envs/roary.yml'
    resources: mem_mb = 32000
    threads: 4
    params: outdir = lambda wildcards: "outputs/roary/" + wildcards.acc + "/"
    shell:'''
    roary -e -n -f {params.outdir} -p {threads} -z {input}
    '''

rule roary_translate:
    input: 'outputs/roary/{acc}/pan_genome_reference.fa' 
    output: 'outputs/roary/{acc}/pan_genome_reference.faa' 
    conda: 'envs/emboss.yml'
    resources:
        mem_mb = 16000
    threads: 2
    shell:'''
    transeq {input} {output}
    '''
