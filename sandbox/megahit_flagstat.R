flagstat <- read_tsv("sandbox/test_megahit_diginorm_nocat/sandbox/multiqc_cdhit_samtools_flagstat.txt")
colnames(flagstat)

ggplot(flagstat, aes(x = mapped_passed_pct)) +
  geom_density() + 
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  labs(x = "percent reads aligned", 
       caption = paste0("mean = ", round(mean(flagstat$mapped_passed_pct, na.rm = T), digits = 1), "\n",
                        "sd = ", round(sd(flagstat$mapped_passed_pct, na.rm = T), digits = 1)))

hist(flagstat$mapped_passed_pct)
flagstat$library <- gsub('_.*', '' , flagstat$Sample)
flagstat$genome <- gsub("^[^_]*_", "", flagstat$Sample)
flagstat$genome <- gsub("^TR_", "", flagstat$genome)

# what was the per genome mapping rate per library?
lib <- flagstat %>%
  group_by(library) %>%
  summarize(mean(mapped_passed_pct))

genome <- flagstat %>%
  group_by(genome) %>%
  summarize(mean(mapped_passed_pct))

ggplot(flagstat %>% filter(genome == "ERS396297_11.fna"), aes(x = mapped_passed_pct)) +
  geom_density() + 
  theme_minimal()

# ERS396297_11 checkm

checkm <- read_tsv("~/Downloads/completeness.tsv", 
                   col_names = c("id", "lineage", "n_genomes", "n_markers", 
                                 "n_marker_sets", "0", "1", "2", "3", "4", "5+", 
                                 "completeness", "contamination", "strain_het"),
                   skip = 1)

ggplot(checkm, aes(x = completeness)) +
  geom_histogram(bins = 50) +
  theme_minimal()

qual <- checkm %>%
  filter(completeness > 80) %>%
  filter(strain_het < 10) %>%
  filter(contamination < 10)

qual_good_names <- gsub("_.*", "", qual$id)
qual_good_names

tmp <- flagstat %>% 
  filter(library %in% qual_good_names) %>%
  filter(genome == "ERS396297_11.fna")
View(tmp)

genomes_to_check <- tmp$Sample
cat(genomes_to_check)

flagstat_cdhit <- read_tsv("~/Downloads/multiqc_samtools_flagstat_cdhit.txt")
hist(flagstat_cdhit$mapped_passed_pct)
flagstat_cdhit$library <- gsub('_.*', '' , flagstat_cdhit$Sample)
flagstat_cdhit$genome <- gsub("^[^_]*_", "", flagstat_cdhit$Sample)
flagstat_cdhit$genome <- gsub("^TR_", "", flagstat_cdhit$genome)


# combine and assess
flagstat <- flagstat %>% 
  select(library, genome, Sample, mapped_passed, flagstat_total, mapped_passed_pct)

flagstat_cdhit <- flagstat_cdhit %>% 
  select(library, genome, Sample, mapped_passed, flagstat_total, mapped_passed_pct)

all <- full_join(flagstat, flagstat_cdhit, by = c("library", "genome", "Sample"))

all$delta_pct <- all$mapped_passed_pct.y - all$mapped_passed_pct.x
mean(all$delta_pct, na.rm = T)



by_genome <- all %>% 
  filter(!is.na(delta_pct)) %>%
  filter(!is.na(mapped_passed_pct.y)) %>%
  select(genome, delta_pct, mapped_passed_pct.x, mapped_passed_pct.y) %>%
  group_by(genome) %>%
  summarise_all(list(mean))
