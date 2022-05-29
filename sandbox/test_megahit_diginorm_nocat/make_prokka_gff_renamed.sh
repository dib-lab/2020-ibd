cd ~/github/2020-ibd/sandbox/test_megahit_diginorm_nocat/prokka
for indir in *
do
bn=$(basename $indir)
cp ${indir}/GCF_900036035.1_RGNV35913_genomic.fna.gff ../prokka_gff_renamed/${bn}_GCF_900036035.1_RGNV35913_genomic.fna.gff
done
