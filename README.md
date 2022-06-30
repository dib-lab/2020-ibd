# ibd
Analysis of publicly available metagenomic sequencing data from humans with IBD.

## Purpose

Implements analyses written up here: https://dib-lab.github.io/2021-paper-ibd/

## Getting started with this repository

The bulk of this pipeline is executed with the Snakefile.

This repository uses conda to manage software installations. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

```
conda env create --name ibd --file environment.yml
conda activate ibd
```

An example command to run on cluster with slurm is included below (dry run):
```
snakemake -j 33 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=400000 --cluster "sbatch -t 1440 -J ibdsmk -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```


## Improvements I would make if I were to do it all again

A list, in no particular order.
Some are trivial, some are a big deal

+ Make the preprocessing steps host remove with bbduk -> quality trim with fastp -> kmer trim with khmer (see https://github.com/dib-lab/distillerycats/blob/main/distillerycats/Snakefile)
+ Update from sourmash compute to sourmash sketch
+ Make spacegraphcats rules proper snakemake checkpoints
+ Upgrade to GTDB rs207 (or whatever GTDB is most current)
+ Do all interpretation analysis in jupyter notebooks (even for R code) so figures are inline with code and easier to keep track of (note I converted most analyses over to notebooks).

