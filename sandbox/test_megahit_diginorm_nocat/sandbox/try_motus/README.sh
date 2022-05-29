

for infolder in ../../../../outputs/sgc_genome_queries/*_k31_r1_search_oh0/ 
do
    bn=$(basename ${infolder} _k31_r1_search_oh0)
    echo "Running motus on ${bn}"
    motus map_snv -s ${infolder}/GCF_900036035.1_RGNV35913_genomic.fna.gz.cdbg_ids.reads.fa.gz -t 8 > map_snv_out/${bn}.bam
done

motus snv_call -d map_snv_out -o snv_call_out 

motus snv_call -fb 40 -fd 2.0 -fp .20 -fc 2 -d map_snv_out -o snv_call_out2

motus snv_call -fb 10 -fd 1 -fp .1 -fc 1 -d map_snv_out -o snv_call_out3
