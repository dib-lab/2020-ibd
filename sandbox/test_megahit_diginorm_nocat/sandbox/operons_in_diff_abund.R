library(dplyr)
library(purrr)
library(readr)

info <- read_tsv("inputs/working_metadata.tsv") %>%
  select(library_name, study_accession, diagnosis) %>%
  distinct()

eggnog_files <- list.files("sandbox/test_megahit_diginorm_nocat/eggnog", 
                           ".annotations$", full.names = T)

eggnog <- eggnog_files %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, .id = "genome", skip = 3) %>%
  mutate(genome = gsub(".emapper.annotations", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/eggnog/", "", genome)) %>%
  rename(query_name = "#query_name")

corncob_files <- list.files("sandbox/test_megahit_diginorm_nocat/corncob",
                            ".sig_ccs.tsv$", full.names = T) 
corncob <- corncob_files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "genome") %>%
  mutate(genome = gsub("_sig_ccs.tsv", "", genome)) %>%
  mutate(genome = gsub("sandbox\\/test_megahit_diginorm_nocat\\/corncob/", "", genome))

corncob$orf_num <- gsub(".*_", "", corncob$aa_seq)
corncob$library_name <- gsub("_+[^_]+$", "", corncob$aa_seq)

# corncob_orf <- corncob %>%
#   filter(mu == "mu.diagnosisCD") %>%
#   filter(estimate < 0) %>%
#   filter(genome == "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna") %>%
#   arrange(library_name, orf_num) %>%
#   mutate(orf_num = as.numeric(orf_num)) %>%
#   mutate(orf_diff = orf_num - lag(orf_num))
#   #filter(orf_diff == 1)
# 
# corncob_orf$orf_diff <- ifelse(is.na(corncob_orf$orf_diff), 0, corncob_orf$orf_diff)
# corncob_orf$orf_stretch <- sequence(rle(as.character(corncob_orf$orf_diff))$lengths)
# 
# tmp <- corncob_orf %>%
#   filter(orf_diff == 1) 
# View(tmp)
# dim(tmp)
# 
# View(corncob_orf)




get_orf_info_cd_down <- function(corncob = corncob, goa){
  corncob_orf <- corncob %>%
    filter(mu == "mu.diagnosisCD") %>%
    filter(estimate < 0) %>%
    filter(genome == goa) %>%
    arrange(library_name, orf_num) %>%
    mutate(orf_num = as.numeric(orf_num)) %>%
    mutate(orf_diff = orf_num - lag(orf_num))
  
  corncob_orf$orf_diff <- ifelse(is.na(corncob_orf$orf_diff), 0, corncob_orf$orf_diff)
  corncob_orf$orf_stretch <- sequence(rle(as.character(corncob_orf$orf_diff))$lengths)
  
  tmp <- corncob_orf %>%
    filter(orf_diff == 1) 
  print(paste0("Dimensions of ORF df: ", dim(tmp)))
  print(paste0("Max consecutive ORFs: ", max(tmp$orf_stretch)))
}
get_orf_info_cd_up <- function(corncob = corncob, goa){
  corncob_orf <- corncob %>%
    filter(mu == "mu.diagnosisCD") %>%
    filter(estimate > 0) %>%
    filter(genome == goa) %>%
    arrange(library_name, orf_num) %>%
    mutate(orf_num = as.numeric(orf_num)) %>%
    mutate(orf_diff = orf_num - lag(orf_num))
  
  corncob_orf$orf_diff <- ifelse(is.na(corncob_orf$orf_diff), 0, corncob_orf$orf_diff)
  corncob_orf$orf_stretch <- sequence(rle(as.character(corncob_orf$orf_diff))$lengths)
  
  tmp <- corncob_orf %>%
    filter(orf_diff == 1) 
  print(paste0("Dimensions of ORF df: ", dim(tmp)))
  print(paste0("Max consecutive ORFs: ", max(tmp$orf_stretch)))
  return(tmp)
}
get_orf_info_uc_down <- function(corncob = corncob, goa){
  corncob_orf <- corncob %>%
    filter(mu == "mu.diagnosisUC") %>%
    filter(estimate < 0) %>%
    filter(genome == goa) %>%
    arrange(library_name, orf_num) %>%
    mutate(orf_num = as.numeric(orf_num)) %>%
    mutate(orf_diff = orf_num - lag(orf_num))
  
  corncob_orf$orf_diff <- ifelse(is.na(corncob_orf$orf_diff), 0, corncob_orf$orf_diff)
  corncob_orf$orf_stretch <- sequence(rle(as.character(corncob_orf$orf_diff))$lengths)
  
  tmp <- corncob_orf %>%
    filter(orf_diff == 1) 
  print(paste0("Dimensions of ORF df: ", dim(tmp)))
  print(paste0("Max consecutive ORFs: ", max(tmp$orf_stretch)))
}
get_orf_info_uc_up <- function(corncob = corncob, goa){
  corncob_orf <- corncob %>%
    filter(mu == "mu.diagnosisUC") %>%
    filter(estimate > 0) %>%
    filter(genome == goa) %>%
    arrange(library_name, orf_num) %>%
    mutate(orf_num = as.numeric(orf_num)) %>%
    mutate(orf_diff = orf_num - lag(orf_num))
  
  corncob_orf$orf_diff <- ifelse(is.na(corncob_orf$orf_diff), 0, corncob_orf$orf_diff)
  corncob_orf$orf_stretch <- sequence(rle(as.character(corncob_orf$orf_diff))$lengths)
  
  tmp <- corncob_orf %>%
    filter(orf_diff == 1) 
  print(paste0("Dimensions of ORF df: ", dim(tmp)))
  print(paste0("Max consecutive ORFs: ", max(tmp$orf_stretch)))
  return(tmp)
}

for(gnome in unique(corncob$genome)){
  print(gnome)
  get_orf_info_cd_up(corncob = corncob, goa = gnome)
}

## No up UC
## Down UC-------
# [1] "SRS075078_49.fna"
# [1] "Dimensions of ORF df: 686" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 12"
# [1] "SRR5127401_bin.3.fa"
# [1] "Dimensions of ORF df: 677" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 16"
# [1] "ERS608524_37.fna"
# [1] "Dimensions of ORF df: 646" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 11"
# [1] "ERS537218_9.fna"
# [1] "Dimensions of ORF df: 126" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 5"
# [1] "ERS235603_16.fna"
# [1] "Dimensions of ORF df: 116" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 5"

## Up CD-----------
# [1] "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna"
# [1] "Dimensions of ORF df: 2186" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 17"
# 
# [1] "GCF_900036035.1_RGNV35913_genomic.fna"
# [1] "Dimensions of ORF df: 717" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 11"

cbolt_orf_cd_up <- get_orf_info_cd_up(corncob = corncob, goa =  "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna")
View(cbolt_orf_cd_up)
rgnav_orf_cd_up <- get_orf_info_cd_up(corncob = corncob, goa =  "GCF_900036035.1_RGNV35913_genomic.fna")
View(rgnav_orf_cd_up)
## Down CD-----------
# [1] "ERS235530_10.fna"
# [1] "Dimensions of ORF df: 889" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 9"
# [1] "ERS235531_43.fna"
# [1] "Dimensions of ORF df: 1520" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 8"
# [1] "ERS235603_16.fna"
# [1] "Dimensions of ORF df: 870" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 9"
# [1] "ERS396297_11.fna"
# [1] "Dimensions of ORF df: 729" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 6"
# [1] "ERS396519_11.fna"
# [1] "Dimensions of ORF df: 737" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 9"
# [1] "ERS473255_26.fna"
# [1] "Dimensions of ORF df: 1301" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 7"
# [1] "ERS537218_9.fna"
# [1] "Dimensions of ORF df: 1345" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 10"
# [1] "ERS537328_30.fna"
# [1] "Dimensions of ORF df: 1682" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 6"
# [1] "ERS537353_12.fna"
# [1] "Dimensions of ORF df: 372" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 5"
# [1] "ERS608524_37.fna"
# [1] "Dimensions of ORF df: 1852" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 14"
# [1] "ERS608576_22.fna"
# [1] "Dimensions of ORF df: 919" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 16"
# [1] "GCF_000371685.1_Clos_bolt_90B3_V1_genomic.fna"
# [1] "Dimensions of ORF df: 72" "Dimensions of ORF df: 13"
# [1] "Max consecutive ORFs: 3"
# [1] "GCF_000508885.1_ASM50888v1_genomic.fna"
# [1] "Dimensions of ORF df: 213" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 5"
# [1] "GCF_001405615.1_13414_6_47_genomic.fna"
# [1] "Dimensions of ORF df: 1122" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 21"
# [1] "GCF_900036035.1_RGNV35913_genomic.fna"
# [1] "Dimensions of ORF df: 215" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 6"
# [1] "LeChatelierE_2013__MH0074__bin.19.fa"
# [1] "Dimensions of ORF df: 1037" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 9"
# [1] "LiJ_2014__O2.UC28-1__bin.61.fa"
# [1] "Dimensions of ORF df: 753" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 7"
# [1] "LiSS_2016__FAT_DON_8-22-0-0__bin.28.fa"
# [1] "Dimensions of ORF df: 1778" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 9"
# [1] "QinJ_2012__CON-091__bin.20.fa"
# [1] "Dimensions of ORF df: 1372" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 7"
# [1] "SRR4305229_bin.5.fa"
# [1] "Dimensions of ORF df: 1096" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 10"
# [1] "SRR5127401_bin.3.fa"
# [1] "Dimensions of ORF df: 1335" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 22"
# [1] "SRR5558047_bin.10.fa"
# [1] "Dimensions of ORF df: 1059" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 13"
# [1] "SRR6028281_bin.3.fa"
# [1] "Dimensions of ORF df: 678" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 10"
# [1] "SRS075078_49.fna"
# [1] "Dimensions of ORF df: 1573" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 23"
# [1] "SRS103987_37.fna"
# [1] "Dimensions of ORF df: 1408" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 14"
# [1] "SRS104400_110.fna"
# [1] "Dimensions of ORF df: 1219" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 9"
# [1] "SRS143598_15.fna"
# [1] "Dimensions of ORF df: 1263" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 16"
# [1] "SRS1719112_8.fna"
# [1] "Dimensions of ORF df: 324" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 7"
# [1] "SRS1719498_9.fna"
# [1] "Dimensions of ORF df: 1477" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 11"
# [1] "SRS1735645_19.fna"
# [1] "Dimensions of ORF df: 1237" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 11"
# [1] "SRS476209_42.fna"
# [1] "Dimensions of ORF df: 1152" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 20"
# [1] "VatanenT_2016__G80445__bin.9.fa"
# [1] "Dimensions of ORF df: 1298" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 7"
# [1] "XieH_2016__YSZC12003_37172__bin.63.fa"
# [1] "Dimensions of ORF df: 240" "Dimensions of ORF df: 13" 
# [1] "Max consecutive ORFs: 7"
# [1] "ZeeviD_2015__PNP_Main_232__bin.27.fa"
# [1] "Dimensions of ORF df: 1020" "Dimensions of ORF df: 13"  
# [1] "Max consecutive ORFs: 13"