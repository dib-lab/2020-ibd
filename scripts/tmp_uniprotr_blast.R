library(readxl)
library(UniprotR)

destfile <- "inputs/pnas.1904099116.sd02.xlsx"
url <- "https://www.pnas.org/highwire/filestream/867579/field_highwire_adjunct_files/2/pnas.1904099116.sd02.xlsx"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto") 
}
supsd2 <- readxl::read_xlsx(destfile)
uniprot <- supsd2$`Uniprot accession`[!is.na(supsd2$`Uniprot accession`)]

# will automatically append .fasta to the end of charcter specific in FileName
GETSeqFastaUniprot(uniprot, FilePath = ".", FileName = "tmp")


# run blast ---------------------------------------------------------------

# makeblastdb -in tmp.fasta -dbtype prot
# blastx -query ../../outputs/sgc_pangenome_catlases_corncob_sequences/GCF_008121495.1_CD_increased_contigs.fa -db tmp.fa.fasta -out out.tsv -outfmt 6


# analyze blast results ---------------------------------------------------

blast <- read_tsv("sandbox/rgnv_blast/out.tsv", 
                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

tmp <- blast %>%
  filter(length > 33) %>%
  filter(pident > 80)
View(tmp)

dom_info <- read_tsv("outputs/sgc_pangenome_catlases/GCF_008121495.1_k31_r10_abund/dom_info.tsv")
cdbg_to_pieces <- read_csv("outputs/sgc_pangenome_catlases/GCF_008121495.1_k31_r10/cdbg_to_pieces.csv")
cdbg_to_pieces_filt <- cdbg_to_pieces %>%
  filter(cdbg_node %in% blast$qseqid)

tmp <- blast %>%
  filter(pident >= 100) %>%
  filter(length >= 10) %>%
  left_join(cdbg_to_pieces_filt, by = c("qseqid" = "cdbg_node")) %>%
  group_by(sseqid, dominator) %>%
  tally()
View(tmp)

length(unique(tmp$dominator))