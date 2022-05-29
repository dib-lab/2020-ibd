# This folder contains the original sgc extract_reads R. gnavus results.
# These reads are "wrong" bc these files were generated when the catlas
# still had a bug in it. However, these files are still useful for testing
# purposes, so they were saved here before all of the sgc output files
# were deleted.

for indir in ~/github/2020-ibd/outputs/sgc_genome_queries/*_r1_search_oh0
do
bn=$(basename $indir _k31_r1_search_oh0)
cp ${indir}/GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz ${bn}_GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz
done 

# make sigs for comp
for infile in *fa.gz
do 
bn=$(basename ${infile} .fa.gz)
sourmash compute -k 31 --scaled 2000 -o ${bn}_k31_scaled2000.sig ${infile}
done
mkdir sourmash
mv *sig sourmash

sourmash compare -k 31 -o rgnv_k31_scaled2000.comp --csv rgnv_k31_scaled2000.csv *scaled2000.sig

