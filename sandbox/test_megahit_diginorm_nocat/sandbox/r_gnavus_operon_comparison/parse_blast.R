library(readr)
library(dplyr)
library(ggplot2)

blast <- read_csv("sandbox/test_megahit_diginorm_nocat/sandbox/r_gnavus_operon_comparison/PA98Y4VC114-Alignment-HitTable.csv",
                  col_names = c("qseqid", "sseqid", "pident", "length",
                                "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                                "evalue", "bitscore"))
blast <-blast %>%
  group_by(qseqid) %>%
  arrange(desc(bitscore)) %>%
  top_n(n = 1) %>%
  filter(pident > 60)

write.table(blast$sseqid, "sandbox/test_megahit_diginorm_nocat/sandbox/r_gnavus_operon_comparison/best_match_seqnames.txt",
            col.names = F, row.names = F, quote = F)
