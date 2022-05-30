for infile in *sig
do 
    j=$(basename $infile .scaled2k.sig)
    
    sourmash signature intersect -o ${j}_filt.sig -A ${infile} -k 31 ${infile} ../../01_greater_than_one_filt_sigs/greater_than_one_count_hashes.sig
    
    sourmash signature rename -o ${j}_filt_named.sig -k 31 ${j}_filt.sig ${j}_filt
done