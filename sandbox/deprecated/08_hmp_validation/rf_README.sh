
## Link in raw signatures
mkdir hmp_validation
cd hmp_validation
cd raw_sigs

ln -s ~/github/cosmo-kmers/outputs/sigs/*sig .
cd ..

## generate a text file of the ~14,000 hashes retained after rf vita variable
## selection: make_rf_vita_hashes.R. Produces rf_vita_hashes.txt. 

conda activate sourmash # sourmash version 2.3.0

## filter raw signatures to ~14,000 hashes from rf vita variable selection.
python ../../scripts/hashvals-to-signature.py -o rf_vita_hashes.sig -k 31 --scaled 2000 --name rf_vita_hashes_k31 --filename rf_vita_hashes.txt rf_vita_hashes.txt

mkdir filt_sigs/
for infile in raw_sigs/*sig
do
j=$(basename $infile .scaled2k.sig)
sourmash signature intersect -o filt_sigs/${j}_rf_vita_filt.sig -A ${infile} -k 31 ${infile} rf_vita_hashes.sig
done

mkdir filt_sig_csvs
for infile in filt_sigs/*sig
do
j=$(basename $infile .sig)
python ../../scripts/sig_to_csv.py ${infile} filt_sig_csvs/${j}.csv
done
