cd blast_against_minot_willis
11187  ln -s ../../cdhit/GCF_900036035.1_RGNV35913_genomic.fna.cdhit.faa .
11189  blastp -query GCF_900036035.1_RGNV35913_genomic.fna.cdhit.faa -subject IBD-associated-genes.fastp -outfmt 6 -out rgvn-ibd.tab

ln -s ../../cdhit/GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna.cdhit.faa .
blastp -query GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna.cdhit.faa -subject IBD-associated-genes.fastp -outfmt 6 -out cbolt-ibd.tab
