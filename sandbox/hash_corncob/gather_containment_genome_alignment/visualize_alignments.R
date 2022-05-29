# https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/

library(ggplot2)
library(dplyr)
library(scales)
library(plotly)
library(magrittr)
library(readr)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(seqinr)

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

# cbolt -------------------------------------------------------------------

coords1 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/c_bolt/GCF_000371525.1_Clos_clos_90A7_VS_GCF_000154365.1_ASM15436v1_genomic.fna_filter.delta")
coords1
coords1$query <- "ASM15436v1"
coords2 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/c_bolt/GCF_000371525.1_Clos_clos_90A7_VS_GCF_000371645.1_Clos_bolt_90B8_V1_genomic.fna_filter.delta")
coords2$query <- "Clos_bolt_90B8"
coords3 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/c_bolt/GCF_000371525.1_Clos_clos_90A7_VS_GCF_000371665.1_Clos_bolt_90B7_V1_genomic.fna_filter.delta")
coords3$query <- "Clos_bolt_90B7"
coords4 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/c_bolt/GCF_000371525.1_Clos_clos_90A7_VS_GCF_000371725.1_Clos_bolt_90A5_V1_genomic.fna_filter.delta")
coords4$query <- "Clos_bolt_90A5"
coords5 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/c_bolt/GCF_000371525.1_Clos_clos_90A7_VS_GCF_001078425.1_Clos_bolt_WAL_14578_V1_genomic.fna_filter.delta")
coords5$query <- "Clos_bolt_WAL_14578"
coords <- rbind(coords1, coords2, coords3, coords4, coords5)

coords %<>% mutate(similarity=1-error/abs(qe-qs))
#mumgp.filt %<>% mutate(similarity=1-error/abs(qe-qs))

p <- ggplot(coords %>% filter(rid == "NZ_KB851045.1") , aes(x=rs, xend=re, y=similarity, yend=similarity, linetype = query, color = query)) + 
  geom_segment(alpha = .9) +
  facet_wrap(~query, ncol = 1) +
  theme_bw() + 
  labs(x = 'reference sequence', y= 'similarity') +
  ylim(.95, 1) +
  scale_x_continuous(labels = comma)
p
ggplotly(p)

# rgnv --------------------------------------------------------------------


coords1 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_BackhedF_2015__SID258_4M__bin.6.fa_filter.delta")
coords1$query <- "backhed"
coords2 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_ChengpingW_2017__AS130raw__bin.7.fa_filter.delta")
coords2$query <- "chengping"
coords3 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_ERR1293630_bin.1.fa_filter.delta")
coords3$query <- "bin1"
coords4 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_ERR525797_bin.4.fa_filter.delta")
coords4$query <- "bin4"
coords5 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_GCA_900036035.1_RGNV35913_genomic.fna_filter.delta")
coords5$query <- "canon"
coords6 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_GCF_009831375.1_ASM983137v1_genomic.fna_filter.delta")
coords6$query <- "asm983"
coords7 <- readDelta("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_VS_SRR3131756_bin.8.fa_filter.delta")
coords7$query <- "bin8"

coords <- rbind(coords1, coords2, coords3, coords4, coords5, coords6, coords7)

coords %<>% mutate(similarity=1-error/abs(qe-qs))

coords <- coords %>%
  mutate(rlen = abs(re - rs)) %>%
  mutate(qlen = abs(qe - qs))

# calculate how many basepairs is shared between the "reference" and each of the strains
coords %>%
  group_by(query) %>%
  summarise(total_shared_bp = sum(qlen)) 


# filter to contig (rid) where all seven of the genomes have overlap
rid_filt <- coords %>%
  select(rid, query) %>%
  distinct() %>%
  group_by(rid) %>%
  tally() %>% 
  filter(n == 7)

rid_filt_strict <- coords %>%
  filter(similarity > .99) %>%
  select(rid, query) %>%
  distinct() %>%
  group_by(rid) %>%
  tally() %>% 
  filter(n == 7)

# calculate how many basepairs is shared by all 7 strains
coords %>% 
  filter(rid %in% rid_filt$rid) %>%
  group_by(query) %>%
  summarise(shared_bp_all = sum(qlen))

p <- ggplot(coords %>%
              filter(rid %in% rid_filt_strict$rid), aes(x=rs, xend=re, y=similarity, yend=similarity, linetype = query, color = query)) + 
  geom_segment(alpha = .9) +
  facet_wrap(~rid, scales = "free_x") +
  theme_bw() + 
  labs(x = 'reference sequence', y= 'similarity') +
  #ylim(.95, 1) +
  scale_x_continuous(labels = comma)
pdf("tmp_rgnv_free_x_99.pdf", height = 40, width = 40)
p
dev.off()

# futher filter to contigs that contain hashes that are important/up in CD
contain_err <- read_csv("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/ERR1293630_bin.1_vs_up_cd_hashes_containment.txt") %>%
  filter(similarity > 0)

p <- ggplot(coords %>%
              filter(rid %in% rid_filt$rid) %>%
              filter(rid %in% contain_err$name), aes(x=rs, xend=re, y=similarity, yend=similarity, linetype = query, color = query)) + 
  geom_segment(alpha = .9) +
  facet_wrap(~rid, scales = "free_x") +
  theme_bw() + 
  labs(x = 'reference sequence', y= 'similarity') +
  ylim(.95, 1) +
  scale_x_continuous(labels = comma)
pdf("tmp_rgnv_free_x.pdf", height = 40, width = 40)
p
dev.off()
ggplotly(p)


# try intersect -----------------------------------------------------------

gr1 <- makeGRangesFromDataFrame(coords1, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr2 <- makeGRangesFromDataFrame(coords2, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr3 <- makeGRangesFromDataFrame(coords3, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr4 <- makeGRangesFromDataFrame(coords4, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr5 <- makeGRangesFromDataFrame(coords5, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr6 <- makeGRangesFromDataFrame(coords6, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr7 <- makeGRangesFromDataFrame(coords7, keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
coord_int <- GenomicRanges::intersect(gr1, gr2)
coord_int <-  GenomicRanges::intersect(coord_int, gr3)
coord_int <-  GenomicRanges::intersect(coord_int, gr4)
coord_int <-  GenomicRanges::intersect(coord_int, gr5)
coord_int <-  GenomicRanges::intersect(coord_int, gr6)
coord_int <-  GenomicRanges::intersect(coord_int, gr7)

df <- GenomicRanges::as.data.frame(coord_int)
View(df)
sum(df$width)


# greater than 98% similarity ---------------------------------------------

gr1 <- makeGRangesFromDataFrame(coords1 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr2 <- makeGRangesFromDataFrame(coords2 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr3 <- makeGRangesFromDataFrame(coords3 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr4 <- makeGRangesFromDataFrame(coords4 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr5 <- makeGRangesFromDataFrame(coords5 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr6 <- makeGRangesFromDataFrame(coords6 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
gr7 <- makeGRangesFromDataFrame(coords7 %>%
                                  mutate(similarity=1-error/abs(qe-qs)) %>%
                                  filter(similarity > .98), keep.extra.columns=FALSE, ignore.strand=T,
                                seqinfo=NULL, seqnames.field=c("rid"), start.field="rs",
                                end.field="re", starts.in.df.are.0based=FALSE)
coord_int <- GenomicRanges::intersect(gr1, gr2)
coord_int <-  GenomicRanges::intersect(coord_int, gr3)
coord_int <-  GenomicRanges::intersect(coord_int, gr4)
coord_int <-  GenomicRanges::intersect(coord_int, gr5)
coord_int <-  GenomicRanges::intersect(coord_int, gr6)
coord_int <-  GenomicRanges::intersect(coord_int, gr7)

df <- GenomicRanges::as.data.frame(coord_int)
# sum(df$width)
# table(df$seqnames %in% contain_err$name)

export.bed(object = coord_int, con ='sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/intersection/intersection_98.bed')

# bedtools getfasta -fi ERR1293630_bin.1.fa -bed intersection_98.bed -name > ERR1293630_bin.1_intersection_98.fa  


winsize = 500
genome = read.fasta(file = "sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/GCF_009831375.1_ASM983137v1_genomic.fna")

for (sequence in genome) {
  
  starts = seq(1, length(sequence) - winsize, by = winsize)
  n = length(starts)
  chunkGCs = numeric(n)
  
  for (i in 1:n) {
    chunk = sequence[starts[i]:(starts[i]+winsize)]
    chunkGC = GC(chunk)
    chunkGCs[i] = chunkGC
  }
  gc_df <- data.frame(window_start = starts, gc = chunkGCs)
  ggplot(gc_df, aes(x = starts, y = gc)) +
    geom_point(size = .5, alpha = .5)+
    theme_minimal()
  
  pdf(file = paste('tmp_gc_plot-', attr(sequence, 'name'), '.pdf', sep = ''))
  plot(starts, chunkGCs, type='b')
  dev.off()
  print(paste(attr(sequence, 'name'), ' done'))
}

read_blast_tab <- function(path){
  blast <- read_tsv(path,
                    col_names = c("qseqid", "sseqid", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                                  "evalue", "bitscore"))
}
# blastn -subject GCF_009831375.1_ASM983137v1_genomic.fna -query ERR1293630_bin.1_intersection_98.fa -out GCF_009831375.1_ASM983137v1_VS_ERR1293630_bin.1_intersection_98_blast.tab -outfmt 6
# prokka ERR1293630_bin.1_intersection_98.fa --outdir ERR1293630_bin.1_intersection_98_prokka --prefix ERR1293630_bin.1_intersection_98 --metagenome --force --locustag ERR1293630_bin.1_int98 --centre X --compliant
# blastn -subject GCF_009831375.1_ASM983137v1_genomic.fna -query ERR1293630_bin.1_intersection_98_prokka/ERR1293630_bin.1_intersection_98.ffn -out GCF_009831375.1_ASM983137v1_VS_ERR1293630_bin.1_intersection_98_genes_blast.tab -outfmt 6
blast <- read_blast_tab("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/intersection/GCF_009831375.1_ASM983137v1_VS_ERR1293630_bin.1_intersection_98_blast.tab")

blast <- blast %>%
  filter(pident > 98) 
 
# DEMARCATE CONTIGS THAT CONTAINED HASHES
blast <- blast %>%
  mutate(qseqid = gsub("\\.::", "", qseqid)) %>%
  mutate(qseqid = gsub(":.*", "", qseqid)) %>%
  filter(qseqid %in% contain_err$name)

gc_df <- gc_df %>%
  mutate(window_end = window_start + 499) %>%
  select(window_start, window_end, gc) %>%
  mutate(feature = "GC") %>%
  mutate(annot = NA) %>%
  mutate(marker = NA)
  

blast <- blast %>%
  mutate(window_start = sstart) %>%
  mutate(window_end = send) %>%
  mutate(feature = "blast") %>%
  mutate(gc = .2) %>%
  mutate(annot = NA) %>%
  mutate(marker = NA) %>%
  select(window_start, window_end, gc, feature, annot, marker)

all_df <- rbind(gc_df, blast)
all_df$feature <- factor(all_df$feature, levels = c("GC", "blast"))
ggplot(all_df, aes(x = window_start, xend = window_end, y = gc, yend = gc, color = feature)) +
  geom_segment(size = 6) +
  theme_minimal() +
  scale_color_manual(values = c(blast = "grey", GC = "black"))


# annotate with prokka annotations
marker_genes <- c("rpsK", "rpsQ", "rpsB", "rpsH", "rpsP", "rpsD", "rpsO", "rpsT", "rpsS", 
                  "rpsE", "rpsC", "rpsI", "rpsJ", "rpsU", "rpsR", "rpsF", "rpsN", "rpsM", 
                  "rpsL", "rpsG", "rplU", "rplQ", "rplP", "rplS", "rplL", "rplD", "rplB", 
                  "rplC", "rplW", "rplE", "rplV", "rplI", "rplM", "rplA", "rplX", "rplJ", "rplR",
                  "rplN", "rplO", "rplT", "rpoB", "rpoA", "rpoC", "rpoN", "rpoD",  "rpmC",  
                  "rpmJ",  "rpmB",  "rpmA",  "rpmE2", "rpmD", "gyrA", "gyrB", "ispH", "rsgA",
                  "lepA", "yqgF", "nusA", "nusG", "uvrC", "secA", "secY", "secG", "recR", "recA",
                  "rumA", "rumA1", "spoU", "rrmJ", "rrmA", "infB", "alaS", "argS", 
                  "ileS", "ctgA", "coaE", "cysS", "dnaA", "dnaG", "dnaX", "engA", "ffh", 
                  "frr", "ftsY", "gmk", "hisS", "infC", "leuS", "ligA", 
                  "pgk", "pheS", "pheT", "prfA", "pyrG", "rbfA", "rnc", "serS", "tsaD",
                  "uvrB", "ybeY", "ychF", "aspS", "ksgA", "fmt", "tilS", "tsf", "trmD",
                  "tig", "smpB", "truB")

prokka_annot <- read_tsv("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/intersection/ERR1293630_bin.1_intersection_98_prokka/ERR1293630_bin.1_intersection_98.tsv")

prokka_annot <- prokka_annot %>%
  mutate(marker =  ifelse(str_detect(string = product, pattern = "Ribosomal"), "marker", 
                          ifelse(gene %in% marker_genes, "marker", "other"))) 
prokka_blast <- read_blast_tab("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/intersection/GCF_009831375.1_ASM983137v1_VS_ERR1293630_bin.1_intersection_98_genes_blast.tab")
prokka <- left_join(prokka_blast, prokka_annot, by = c("qseqid" = "locus_tag")) %>%
  filter(ftype =="CDS") %>%
  mutate(annot = paste0(gene, ": ", product)) %>%
  mutate(gc = .15) %>%
  mutate(feature = "gene") %>%
  mutate(window_start = sstart) %>%
  mutate(window_end = send) %>%
  select(window_start, window_end, gc, feature, annot, marker) %>%
  distinct()

all_df <- rbind(all_df, prokka)

# tblastn -subject GCF_009831375.1_ASM983137v1_genomic.fna -query rumgna_03512-03534.fasta -out GCF_009831375.1_ASM983137v1_VS_inflamm_poly.tab -outfmt 6
inflamm_poly <- read_blast_tab("sandbox/hash_corncob/gather_containment_genome_alignment/rgnv/intersection/GCF_009831375.1_ASM983137v1_VS_inflamm_poly.tab") %>%
  filter(pident > 90) %>%
  mutate(window_start = sstart) %>%
  mutate(window_end = send) %>%
  mutate(feature = "inflamm_poly") %>%
  mutate(gc = .1) %>%
  mutate(annot = NA) %>%
  mutate(marker = NA) %>%
  select(window_start, window_end, gc, feature, annot, marker)

all_df <- rbind(all_df, inflamm_poly)

pdf("gc_map_int_annot_hash_contigs_marked.pdf", height = 5, width = 50)
ggplot(all_df, aes(x = window_start, xend = window_end, y = gc, yend = gc, 
                   color = marker, label = annot)) +
  geom_segment(size = 8) +
  #geom_text(size = .5, angle = 90, vjust = -.2, color = "black") +
  theme_minimal() +
  ylim(-.2, .7) 
  #scale_color_manual(values = c(blast = "grey", GC = "black", gene = "grey"))
dev.off()

plt <- ggplot(all_df, aes(x = window_start, y = gc, 
                          color = marker, label = annot)) +
  geom_point(size = 1, alpha = .9) +
  theme_minimal() 
ggplotly(plt)
