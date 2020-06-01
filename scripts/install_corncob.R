library(devtools)

devtools::install_github("bryandmartin/corncob")
library(corncob)
file.create(snakemake@output[['corncob']])
