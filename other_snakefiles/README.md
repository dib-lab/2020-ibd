This directory contains analyses that didn't make it into the final paper.
Most of the code was written and executed during the first round of analysis after the first time we ran spacegraphcats genome queries on the metagenomes.
Subsequently, we identified a bug is spacegraphcats.
Because we needed to re-run spacegraphcats, we chose to re-build our models as well, combining many models instead of focusing on just one.
At this point, we pivoted our analysis to dominating set differential abundance analysis and annotation directly on the graph, meaning a lot of the analyses we can done previously were no longer relevant.

Most of the snakefile are not stand alone pipelines. 
Some will run if the Snakefile in the main directory has been executed and all of the output files exist.
However, they do document an okay set of rules to acheive different things that at different points in this research project were useful, so I'm keeping them around to document those steps.

+ `singlem_models.snakefile`: These rules build and assess models built on the output of the singlem tool.
SingleM identifies and estimates the abundance of single copy marker genes in short shotgun metagenome sequencing reads.
These sequences are used to infer the taxonomic composition of the metagenome.
The purpose of these rules was to determine if the FracMinHash random forest classifiers only captured signal from marker genes, or if they also captured signal from other sequences in the metagenomes.
For a write up of these results, see my thesis: https://github.com/taylorreiter/2020-dissertation
Chapter 4: https://github.com/taylorreiter/2020-dissertation/blob/main/thesis/04-ibd.Rmd#L155
Rendered document: https://github.com/taylorreiter/2020-dissertation/releases/download/2020_11_18/thesis.pdf
As noted in the snakefile, the models that were built from the spacegraphcats query neighborhood outputs were never re-run after we re-ran spacegraphcats.

+ `orpheum_sgc_nbhds.snakefile`: These rules translate the spacegraphcats genome query neighborhood reads into amino acid sequences in the correct open reading frame.
These rules were run on the new spacegraphcats outputs. 
Orpheum predicts the open reading frame of short sequencing reads and translates the reads into protein space.
Originally, we had planned to do metapangenome analysis on the genome query neighborhoods (see https://github.com/taylorreiter/2021-paper-metapangenomes/).
However, the goal of this analysis was to determine if sequences that were more abundant in IBD (CD or UC) were exclusive to IBD.
I wanted to do this in protein space because it is more permissive to small evolutionary difference (e.g. third base pair wobble).
This ended up being difficult though, because we ran orpheum on the spacegraphcats neighborhoods from each metagenome, and ran differential abundance analysis on the metapangenome species graph (graph that contained sequences from all metagenomes, for one particular species).
Unlike the metagenome queries, it wasn't as clear which specific reads led to each dominating set piece in the metapangenome species graph, and therefore it was hard to get the reads back that underlied the differentially abundant dominating set pieces.
Without reads, I couldn't run orpheum.
What I ended up doing instead was using nucleotide k-mers from the dominating set pieces, because I didn't need to go back to the reads to get the translation.
This ended up serving my purpose fine, but made the orpheum rules obsolete.
 
+ `mgx_assembly_diff_abund_analysis.snakefile`: This snakefile records rules to do a traditional metagenome assembly analysis on the spacegraphcats genome query neighborhoods.
This type of analysis pipeline is replaced by the dominating set differential abundance analysis method paired with catlas annotation.
The draw back of doing this sort of analysis is that many of the reads will not assemble.
Our original analysis demonstrated that while ~80-90% of reads will map back to the assembly, many of the reads that contain the hashes that were most correlated with IBD did not assemble. 
See here for a preliminary description: https://github.com/taylorreiter/2020-dissertation/blob/main/thesis/04-ibd.Rmd#L259
This snakefile was not re-run on the new spacegraphcats results.

+ `roary.snakefile`: this snakefile was a proof of concept/minimum running example that allowed me to demo an approach for collecting sequences to use for multifasta annoations.
The approach matured and was semi-integrated into the main snakefile, and then was pulled back out and put into `unused_generate_multifasta_queries.snakefile`.

+ `unused_generate_multifasta_queries.snakefile`: The most complete way to do multifasta queries for annotations (as of May 2022).
Multifasta queries report the cdbg ids that contain an annotation.
The annotations are fed to the multifasta queries as a multifasta file. 
Each gene in the multifasta file is compared to the cdbg nodes.
To generate the most complete set of potential genes to annotate with, this snakefile takes a combined reference and de novo approach.
It first uses a referenced-based approach to build a pangenome from genomes of the same species in GTDB rs202. 
This creates a representative set of all known genes that belong to a specific species.
The next step assembles the individual genome query neighborhoods, annotates them with prokka, and dereplicates them by clustering with cdhit. 
These two sets of genes are then combined and used for the multifasta query.
The rules in this snakefile were never completely run, so I'm not sure if they contains bugs or not.
I didn't end up using this approach because the reference-based approach alone was so much more light weight and gave me all of the information I needed to understand the high level biology that occurred in this system.
