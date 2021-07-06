# ibd
Analysis of publicly available metagenomic sequencing data from humans with IBD.

Example run command on cluster with slurm

```
snakemake -j 33 --use-conda --rerun-incomplete --restart-times 1 --latency-wait 15 --until loo_validation_seed --resources mem_mb=400000 --default-resources runtime=1440 --cluster "sbatch -t {resources.runtime} -J ibd_smake -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -n
```
