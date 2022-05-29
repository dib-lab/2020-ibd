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
        import pdb; pdb.set_trace()
        return p
#=========================

class Checkpoint_RoaryPrefetchResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {roary_acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self, acc):
        prefetch_csv = f'outputs/genbank/{acc}.x.genbank.gather.csv'
        assert os.path.exists(prefetch_csv)

        genome_accs = []
        with open(prefetch_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               roary_acc = row['name'].split(' ')[0]
               genome_accs.append(roary_acc)
        print(f'loaded {len(genome_accs)} accessions from {prefetch_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'gather_gtdb_rep_to_shared_assemblies';
        # this will trigger exception until that rule has been run.
        #import pdb; pdb.set_trace()
        checkpoints.roary_filter_prefetch_shared_assemblies_vs_refseq.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(w.acc)

        p = expand(self.pattern, roary_acc=genome_accs, **w)
        return p

rule all:
    input: 
        Checkpoint_GatherResults(Checkpoint_RoaryPrefetchResults("outputs/roary_prokka/{{acc}}/{roary_acc}.gff"))

#=======================================
checkpoint gather_gtdb_rep_to_shared_assemblies:
    output: 
        gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    touch {output}
    '''

rule download_shared_assemblies:
    input: 
        gather_grist = "outputs/genbank/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv",
        conf = "inputs/genome-grist-conf.yml"
    output: "genbank_genomes/{acc}_genomic.fna.gz"
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
    output: 
        filt="outputs/genbank/{acc}.x.genbank.gather.csv",
        final_filenames="outputs/genbank/{acc}_roary_input_filenames.txt"
    conda: "envs/tidy.R"
    script: "scripts/roary_filter_prefetch.R"

 
rule roary_make_conf_files_for_genome_grist:
    input: "outputs/genbank/{acc}.x.genbank.gather.csv"
    output: conf="outputs/roary_genome_grist_conf/{acc}_conf.yml"
    run:
        with open(output.conf, 'wt') as fp:
           print(f"""\
sample:
- {wildcards.acc}
outdir: outputs
""", file=fp)

rule roary_genome_grist:
    input:
        prefetch = "outputs/genbank/{acc}.x.genbank.gather.csv",
        conf="outputs/roary_genome_grist_conf/{acc}_conf.yml"
    output: "genbank_genomes/{acc}/{roary_acc}_genomic.fna.gz"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
    '''


rule roary_prokka:
    output: 
        ffn = 'outputs/roary_prokka/{acc}/{roary_acc}.ffn',
        faa = 'outputs/roary_prokka/{acc}/{roary_acc}.faa',
        gff = "outputs/roary_prokka/{acc}/{roary_acc}.gff"
    conda: 'envs/prokka.yml'
    resources:
        mem_mb = 8000
    threads: 2
    params: 
        outdir = lambda wildards: 'outputs/roary_prokka/' + wildcards.acc + "/",
        prefix = lambda wildcards: wildcards.roary_acc,
        gzip = lambda wildcards: "genbank_genomes/" + wildcards.roary_acc + "_genomic.fna.gz"
    shell:'''
    gunzip {input}
    prokka {params.gzip} --outdir {params.outdir} --prefix {params.prefix} --metagenome --force --locustag {params.prefix} --cpus {threads} --centre X --compliant
    mv {params.prefix}.ffn {output.ffn}
    mv {params.prefix}.faa {output.faa}
    gzip {params.gzip}
    '''

# rule roary
