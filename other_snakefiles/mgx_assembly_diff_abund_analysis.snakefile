# This snakefile records rules to do a traditional metagenome assembly analysis on the spacegraphcats genome query neighborhoods.
# This type of analysis pipeline is replaced by the dominating set differential abundance analysis method paired with catlas annotation.
# The draw back of doing this sort of analysis is that many of the reads will not assemble.
# Our original analysis demonstrated that while ~80-90% of reads will map back to the assembly, many of the reads that contain the hashes that were most correlated with IBD did not assemble. 
# See here for a preliminary description: https://github.com/taylorreiter/2020-dissertation/blob/main/thesis/04-ibd.Rmd#L259
# We eventually re-ran spacegraphcats after a bug was fixed, so the exact number would change but this gives an idea of the magnitude of the assembly problem.

# This is not a standalone pipeline.
# It relies on output files that are encoded in the Snakefile.
# These rules could be dropped in to the snakefile and they should run (or at least only need minor modification before they run)

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

