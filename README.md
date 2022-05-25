# ibd
Analysis of publicly available metagenomic sequencing data from humans with IBD.

Example run command on cluster with slurm

```
snakemake -j 33 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=400000 --cluster "sbatch -t 1440 -J ibdsmk -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
