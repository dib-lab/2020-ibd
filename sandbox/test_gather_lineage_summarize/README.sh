cat ../../outputs/gather/*genbank*csv > combined_genbank_results.csv
wget -O gather-to-tax.py https://github.com/dib-lab/sourmash/files/4736660/gather-to-tax.py.txt
chmod +x gather-to-tax.py
gather-to-tax.py species combined_genbank_results.csv /home/irber/sourmash_databases/outputs/lca/lineages/bacteria_genbank_lineage.csv
cat /home/irber/sourmash_databases/outputs/lca/lineages/*_genbank_lineage.csv > all_genbank_lineages.csv

for infile in ../../outputs/gather/*genbank_seed*.csv
do
bn=$(basename $infile .csv)
./gather-to-tax.py $infile all_genbank_lineages.csv > ${bn}_lineage_summarized.tsv
done

# combine genbank taxonomy files
cat /home/irber/sourmash_databases/outputs/lca/lineages/protozoa_genbank_lineage.csv /home/irber/sourmash_databases/outputs/lca/lineages/fungi_genbank_lineage.csv /home/irber/sourmash_databases/outputs/lca/lineages/viral_genbank_lineage.csv > protozoa_fungi_viral_genbank_lineage.csv
# format tessa's gtdb lineage sheet
# add a blank col to gtdb where taxonid is stored
wc -l /group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv
paste -d, <(cut -d, -f-1 /group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv) <(for i in `seq 1 258407`; do echo; done) <(cut -d, -f2- /group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv) > gtdb-rs202.taxonomy_cut.csv
# remove version from accession ID
sed 's/\.[^,]*//' gtdb-rs202.taxonomy_cut.csv > gtdb-rs202.taxonomy_cut_acc.csv
cat gtdb-rs202.taxonomy_cut_acc.csv protozoa_fungi_viral_genbank_lineage.csv > gtdb-rs202-genbank-protozoa-viral-fungi-lineage.csv

for infile in ../../outputs/gather/*gtdb_seed*.csv
do
bn=$(basename $infile .csv)
./gather-to-tax.py $infile gtdb-rs202-genbank-protozoa-viral-fungi-lineage.csv > ${bn}_lineage_summarized.tsv
done
