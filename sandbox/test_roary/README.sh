cp ../../genbank_genomes/GCF_008121495.1_genomic.fna.gz .
conda create -n roary roary sourmash prokka 
conda activate roary
sourmash compute -k 31 --scaled 2000 -o GCF_008121495.1_genomic.sig GCF_008121495.1_genomic.fna.gz
sourmash prefetch -k 31 -o GCF_008121495.1_genomic_prefetch.csv GCF_008121495.1_genomic.sig /home/irber/sourmash_databases/outputs/sbt/refseq-bacteria-x1e6-k31.sbt.zip

# run R script to filter and format for genome grist
# filter_prefetch_for_roary.R
# save file to outputs/genbank/... for consumption by genome-grist
# outputs/genbank/roary_prefetch.x.genbank.gather.csv
cp GCF_008121495.1_genomic_prefetch_filtered.csv outputs/genbank/roary_prefetch.x.genbank.gather.csv 
conda install pip
python -m pip install genome-grist
genome-grist run roary_genome_grist_conf.yml --until make_sgc_conf --nolock

# run prokka
for infile in genbank_genomes/*fna
do 
bn=$(basename ${infile} _genomic.fna)
prokka ${infile} --outdir outputs/prokka/${bn} --prefix ${bn} --metagenome --force --locustag ${bn} --cpus 1 --centre X --compliant
done

# run roary
#mv outputs/prokka/*/*gff outputs/prokka/gff/
roary -f outputs/roary -p 5 -z outputs/prokka/*/*gff
roary -e -n -f outputs/roary -p 5 -z outputs/prokka/*/*gff

roary -e -n -f outputs/roary_with_megahit_and_isolates -p 8 -z outputs/prokka/*/*gff ../test_megahit_diginorm_nocat/prokka_gff_renamed/*gff
roary -e -n -f outputs/roary_with_megahit -p 8 -z ../test_megahit_diginorm_nocat/prokka_gff_renamed/*gff

# count number of basepairs per genbank genome
for infile in genbank_genomes/*fna
do
grep -v ">" ${infile} | wc | awk '{print $3-$1}' >> genome_bp.txt
done

# count number AA in prokka; multiply by 3 to get num bp
for infile in outputs/prokka/*/*faa
do
grep -v ">" $infile | wc | awk '{print $3-$1}' >> coding_bp.txt
done


