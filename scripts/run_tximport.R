library(dplyr)
library(readr)
library(tximport)


# produce quant file list -------------------------------------------------

files <- unlist(snakemake@input[['quant']]) # read in list of quant files
# files <- list.files("sandbox/try_diffex", "quant.sf$", recursive = T, full.names = T)
files <- files[file.size(files) > 0]        # remove empty files
gather_genome <- snakemake@params[['gather_genome']]
files_root <- gsub("outputs\\/nbhd_reads_salmon\\/", "", files) # derive sample name from file name
gsub_string <- paste0("\\/", gather_genome, "_quant\\/quant\\.sf")
files_root <- gsub(gsub_string, "", files_root) # derive sample name from file name

# read in metadata --------------------------------------------------------

# read in sample metadata
info <- read_tsv(snakemake@input[["info"]]) %>%
  select(study_accession, library_name, diagnosis) %>%
  distinct() 

# read in abundtrim lib size info
libsizes <- read_tsv(snakemake@input[['mqc_fastp']], skip = 1,
                  col_names = c("library_name", "passed_filter", "low_quality", 
                                "too_many_N", "too_short")) %>%
  mutate(library_name = gsub("\\.abundtrim", "", library_name)) %>%
  group_by(library_name) %>%
  mutate(libsize = sum(passed_filter, low_quality, too_many_N, too_short)) %>%
  select(library_name, libsize)

# join metadata
info <- left_join(info, libsizes, by = "library_name")

# import counts -----------------------------------------------------------

# generate a column in the metadata info with the path to the quant.sf salmon file
info_salmon <- info %>%
  filter(library_name %in% files_root)
print(head(info_salmon))
info_salmon <- info_salmon[order(match(info_salmon$library_name, files_root)), ] # match order of info to order of files_root
stopifnot(all.equal(info_salmon$library_name, files_root)) # check that order matches
info_salmon$salmon <- files # add col

counts <- tximport(files = info_salmon$salmon, type = "salmon", txOut = T)
count_info <- as.data.frame(counts$counts)
colnames(count_info) <- info_salmon$library_name

# write full counts to file
write_tsv(count_info, path = snakemake@output[["counts"]])
