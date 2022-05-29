library(readr)
library(dplyr)

read_blast_tab <- function(path){
  blast <- read_tsv(path,
                    col_names = c("qseqid", "sseqid", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore"))
}


# rgnv --------------------------------------------------------------------
blast <- read_blast_tab("sandbox/test_megahit_diginorm_nocat/sandbox/blast_against_minot_willis/rgvn-ibd.tab")

blast <- blast %>%
  filter(pident >= 99) %>%
  filter(length >= 100)
dim(blast)
View(blast)



corncob_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob",
                            ".sig_ccs.tsv$", full.names = T) 
corncob_pan <- corncob_files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_sig_ccs.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome)) 

corncob_rgnv <- corncob_pan %>%
  filter(genome == "GCF_900036035.1_RGNV35913_genomic.fna")
View(corncob_rgnv)

table(blast$qseqid %in% corncob_rgnv$aa_seq)


# cbolt -------------------------------------------------------------------

cbolt_blast <- read_blast_tab("sandbox/test_megahit_diginorm_nocat/sandbox/blast_against_minot_willis/cbolt-ibd.tab")

cbolt_blast <- cbolt_blast %>%
  filter(pident >= 99) %>%
  filter(length >= 100)
dim(cbolt_blast)

corncob_cbolt <- corncob_pan %>%
  filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna")


table(cbolt_blast$qseqid %in% corncob_cbolt$aa_seq)
