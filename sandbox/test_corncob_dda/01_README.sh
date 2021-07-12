conda activate sandbox
conda install ntcard

mkdir ntcard_kmer_count
for infile in ../../outputs/abundtrim/*fq.gz
do
bn=$(basename $infile .fq.gz)
ntcard -k31 -c2000 -t 4 -o ntcard_kmer_count/${bn}.freq $infile &> ntcard_kmer_count/${bn}.fstat
done 
