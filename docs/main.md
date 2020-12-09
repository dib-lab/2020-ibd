---
title: IBD Meta-analysis
author:
    - Taylor Reiter
    - Luiz Irber
    - ...
    - Phillip Brooks
    - Alicia Gingrich
    - C. Titus Brown
date: \today{}
geometry: "left=2.5in,right=.2in,top=1in,bottom=1in"
bibliography: bibliography.bib
header-includes:
    - \usepackage{setspace}
    - \singlespacing
    - \usepackage{lineno}
    - \linenumbers
fontfamily: helvet
---

## Introduction

Metagenomics captures the functional potential of microbial communities through DNA sequencing of genes and organisms. 
Metagenomics has been used to profile many human microbial communities, including those that change in or contribute to disease. 
In particular, human gut microbiomes have been extensively characterized for their potential role in diseases such as obesity [@greenblum2012], type II diabetes [@qin2012], colorectal cancer [@wirbel2019], and inflammatory bowel disease [@lloyd2019; @morgan2012; @hall2017; @franzosa2019]. 
Inflammatory bowel disease (IBD) refers to a spectrum of diseases characterized by chronic inflammation of the intestines and is likely caused by host-mediated inflammatory responses at least in part elicited by microorganisms [@kostic2014]. 
However, no causative or consistent microbial signature has been associated with IBD to date. 

[//]: # need a transition sentence here

Statements about biology, determined once computation is all done

[//]: # In many studies, IBD is associated with lower biodiversity [@lewis2015; @moustafa2018; @qin2010], however one study failed to replicate this finding [@lloyd2019]. Some studies have noted lower Firmicutes and higher Gammaproteobacteria in IBD, however one study noted that disease state is not independently associated with microbiome shifts after correcting for other covariates [@morgan2012]. Another study found that a novel *Ruminococcus gnavus* clade is enriched in IBD and that this enrichment would be masked in 16s rRNA amplicon studies given commensurate shifts in other closely-related species [@hall2017]. Additionally, time-course studies of IBD have shown that inter-individual differences contribute more to microbiome variation than disease status [@lloyd2019]. 

[//]: # Both reactive oxygen species and reactive nitrogen species are likely important in IBD [@keshavarzian2003; @winter2014]. 

Although there is no consistent taxonomic or functional trend in the gut microbiome associated with IBD diagnosis, metagenomic studies conducted unto this point have left substantial portions of data unanalyzed. 
Reference-based pipelines commonly used to analyze metagenomic data from IBD cohorts such as HUMANn2 characterize on average 31%-60% of reads from the human gut microbiome metagenome, as many reads do not closely match sequences in reference databases [@franzosa2014; @lloyd2019]. 
To combat this issue, reference-free approaches like *de novo* assembly and binning are used to generate metagenome-assembled genome bins (MAGs) that represent species-level composites of closely related organisms in a sample.
However, *de novo* approaches fail when there is low-coverage of or high strain variation in gut microbes, or with sequencing error [@olson2017]. 
Even when performed on a massive scale, an average of 12.5% of reads fail to map to all *de novo* assembled organisms from human microbiomes [@pasolli2019].

Here we perform a meta-analysis of six studies of IBD gut metagenome cohorts comprising 260 CD, 132 UC and 213 healthy controls (see **Table {@tbl:cohorts}**) [@lloyd2019; @lewis2015; @hall2017; @franzosa2019; @gevers2014; @qin2010]. 
First, we re-analyzed each study using a consistent k-mer-based, reference-free approach. 
We demonstrate that diagnosis accounts for a small but significant amount of variation between samples. 
Next, we used random forests to predict IBD diagnosis and to determine the k-mers that are predictive of UC and CD. 
Then, we use compact de Bruijn graph queries to reassociate k-mers with sequence context and perform taxonomic and functional characterization of these sequence neighborhoods. 
We find that strain variation is important (ADD MORE HERE AFTER CORNCOB).
Our analysis pipeline is lightweight and is extensible to other association studies in large metagenome sequencing cohorts.

## Results

| **Cohort**  |**Cohort names**          | **Country**      |**Total**|**CD**|**UC**|**nonIBD**|**Reference**|
|-------------|--------------------------|------------------|-------|-----|-----|--------|--------------|
| iHMP        | IBDMDB                   | USA              | 106   | 50  | 30  | 26     | [@lloyd2019] |
| PRJEB2054   | MetaHIT                  | Denmark, Spain   | 124   | 4   | 21  | 99     | [@qin2010]   |
| SRP057027   | NA                       | Canada, USA      | 112   | 87  | 0   | 25     | [@lewis2015] |
| PRJNA385949 | PRISM, STiNKi            | USA              | 17    | 9   | 5   | 3      | [@hall2017]  |
| PRJNA400072 | PRISM, LLDeep, and NLIBD | USA, Netherlands | 218   | 87  | 76  | 55     | [@franzosa2019] |
| PRJNA237362 | RISK                     | North America    | 28    | 23  | 0   | 5      | [@gevers2014] |
| Total       |                          |                  | 605   | 260 | 132 | 213    |               |
Table: Six IBD cohorts used in this meta-analysis. {#tbl:cohorts}

### Annotation-free approach for meta-analysis of IBD metagenomes.

Given that both reference-based and *de novo* methods suffer from substantial and biased loss of information in the analysis of metagenomes [@thomas2019multiple; @breitwieser2019review], we sought a reference- and assembly-free pipeline to fully characterize gut metagenomes of IBD patients (**Figure {@fig:pipeline}**). 
K-mers, words of length *k* in nucleotide sequences, have previously been used for annotation-free characterization of sequencing data (reviewed by @rowe2019levee).
K-mers are suitable for metagenome analysis because they do not need to be present in reference databases to be included in analysis, and because they capture information from reads even when there is low coverage or high strain variation that preclude assembly. 
In particular, scaled MinHash sketching produces compressed representations of k-mers in a metagenome while retaining the sequence diversity in a sample [@pierce2019].
Importantly, this approach creates a consistent set of hashes across samples by retaining the same hashes when the same k-mers are observed. 
This enables comparisons between metagenomes.
Given these attributes, we use scaled MinHash sketches to perform metagenome-wide k-mer association with IBD-subtype. 
We refer to scaled MinHash sketches as *signatures*, and to each subsampled k-mer in a signature as a *hash*. 

We also implemented a consistent preprocessing pipeline to remove erroneous sequences that could falsely deflate similarity between samples. 
We removed adapters, human DNA, and erroneous k-mers, and filtered signatures to retain hashes that were present in multiple signatures. 
These preprocessing steps removed hashes that were likely to be errors while keeping hashes that were real but low abundance. 
7,376,151 hashes remained after preprocessing and filtering. 


![Comparison of common metagenome analysis techniques with the method used in this paper. Metagenomes consist of short (\~50-300 bp) reads derived from sequencing DNA from environmental samples. **A** Reference-based metagenomic analysis. Reads are compared to genomes, genes, or proteins in reference databases to determine the presence and abundance of organisms and proteins in a sample. Unmapped reads are typically discarded from downstream analysis. **B** *De novo* metagenome analysis. Overlapping reads are assembled into longer contiguous seqeunces (\~500bp-150kbp, [@vollmers2017comparing]) and binned into metagenome-assembled genome bins. Bins are analyzed for taxonomy, abundance, and gene content. Reads that fail to assemble and contigs that fail to bin are usually discarded from downstream analysis. **C** Annotation-free approach for meta-analysis of metagenomes. We decompose reads into k-mers and subsample these k-mers, selecting k-mers that evenly represent the sequence diversity within a sample. We then identify interesting k-mers using random forests, and associate these k-mers with genomes in reference databases. Meanwhile, we construct a compact de Bruijn graph (cDBG) that contains all k-mers from a metagenome. We query this graph with known genomes that contain our interesting k-mers to recover sequence diversity nearby our query sequences in the cDBG. In the colored cDBG, light grey nodes indicate nodes that contain at least one identical k-mer to the query, while nodes outlined in orange indicate the nearby sequences recovered via cDBG queries. The combination of all orange nodes produces a sample-specific pangenome that represents the strain variation of closely-related organisms within a single metagenome. We repeat this process for all metagenomes and generate a single metapangenome depicted in orange, blue, and pink.](figures/fig1.pdf){#fig:pipeline}

### K-mers capture variation due to disease subtype

We first sought to understand whether variation due to IBD diagnosis is detectable in gut metagenomes.  
We calculated pairwise distance matrices using jaccard distance and cosine distance between filtered signatures, where jaccard distance captured sample richness and cosine distance captured sample diversity. 
We performed principle coordinate analysis and PERMANOVA with these distance matrices (**Figure {@fig:comp-plts}**), using the variables study accession, diagnosis, library size, and number of hashes in a filtered signature (**Table {@tbl:permanova}**). 
Number of hashes in a filtered signature accounts for the highest variation, possibly reflecting reduced diversity in stool metagenomes of CD and UC patients (reviewed in [@schirmer2019microbial]). 
Study accounts for the second highest variation, emphasizing that technical artifacts can introduce biases with strong signals.
Diagnosis accounts for a similar amount of variation as study, indicating that there is a small but detectable signal of IBD subtype in stool metagenomes.   

![Principle coordinate analysis of metagenomes from IBD cohorts performed on filtered signatures. **A** Jaccard distance. **B** Angular distance.](figures/fig2.pdf){#fig:comp-plts}
  

| **Variable** | **Jaccard distance** | **Angular distance**|
|--------------|----------------------|---------------------|
| Number of hashes | 9.9%\*           | 6.2%\*|
| Study accession  | 6.6%\*           | 13.5%\*|
| Diagnosis        | 6.2%\*           | 3.3%\* |
| Library size     | 0.009%\*         | 0.01%\*|
Table: Results from PERMANOVA performed on Jaccard and Angular distance matrices. Number of hashes refers to the number of hashes in the filtered signature, while library size refers to the number of raw reads per sample. \* denotes p < .001. {#tbl:permanova}

### Hashes are weakly predictive of IBD subtype

To evaluate whether the variation captured by diagnosis is predictive of IBD disease subtype, we built random forests classifiers to predict CD, UC, or non-IBD.
We used random forests because of the interpretability of feature importance via variable importance measurments.
We used a leave-one-study-out cross-validation approach where we built and optimized a classifier using five cohorts and validated on the sixth.  
Given the high-dimensional structure of this dataset (e.g. many more hashes than samples), we first used the vita method to select predictive hashes in the training set [@janitza2018; @degenhardt2017]. 
Vita variable selection is based on permuation of variable importance, where p-values for variable importance are calculated against a null distribution that is built from variables that are estimated as non-important [@janitza2018].
This approach retains important variables that are correlated [@janitza2018; @seifert2019], which is desirable in omics-settings where correlated features are often involved in a coordinated biological response, e.g. part of the same operon, pathways, or genome [@stuart2003gene; @sabatti2002co]. 
Variable selection reduced the number of hashes used in each model to 29,264-41,701 (**Table {@tbl:varselhashes}**). 
Using this reduced set of hashes, we then optimized each random forests classifier on the training set, producing six optimized models.
We validated each model on the left-out study.
The accuracy on the validation studies ranged from 49.1%-75.9% (**Figure {@fig:acc-plts}**), outperforming a previously published model built on metagenomic data alone [@franzosa2019].

|**Validation study**|**Selected hashes** |**Percent of total hashes**|
|--------------------|--------------------| --------------------------|
| PRJNA385949 | 41701 | 0.57% |
| PRJNA237362 | 40726 | 0.55% |
| iHMP        | 39628 | 0.54% |
| PRJEB2054   | 35343 | 0.48% |
| PRJNA400072 | 32578 | 0.44% |
| SRP057027   | 29264 | 0.40% |
Table: Number of predictive hashes after variable selection for each of 6 classifiers. Classifiers are labelled by the validation study that was held out from training. {#tbl:varselhashes}

We next sought to understand whether there was a consistent biological signal captured among classifiers by evaluating the fraction of shared hashes between models.
We intersected each set of hashes used to build each optimized classifier (**Figure {@fig:acc-plts}**).
Nine hundred thrity two hashes were shared between all classifiers, while 3,859 hashes were shared between at least five studies.
The presence of shared hashes between classifiers indicates that there is a weak but consistent biological signal for IBD subtype between cohorts.      

Shared hashes accounted for 2.8% of all hashes used to build the optimized classifiers.
If shared hashes are predictive of IBD subtype, we would expect that these hashes would account for an outsized proportion of variable importance in the optimized classifiers.
After normalizing variable importance across classifiers, 40.2% of the total variable importance was held by shared between hashes, with 21.5% attributable to the 932 hashes shared between all six classifiers. 
This indicates that shared hashes contribute a large fraction of predictive power for classification of IBD subtype. 

![Random forest classifiers weakly predict IBD subtype. **A** Accuracy of leave-one-study-out random forest classifiers on training and validation sets. The validation study is on the x axis. **B** Confusion matrices depicting performance of each leave-one-study-out random forest classifier on the validation set. **C** Upset plot depicting intersections of sets of hashes as well as the cumulative normalized variable importance of those hashes in the optimized random forest classifiers. Each classifier is labelled by the left-out validation study.](./figures/fig3.pdf){#fig:acc-plts}

Many hashes were identifiable when compared against all microbial genomes in GenBank, as well as metagenome-assembled genomes from three recent *de novo* assembly efforts from human microbiome metagenomes [@pasolli2019; @nayfach2019; @almeida2019]. 
77.7% of hashes from all classifiers were identifiable and anchored to 1,161 genomes (**Figure {@fig:genomes-plt} A**). 
In contrast, 69.4% of shared hashes anchored to 41 genomes (**Figure {@fig:genomes-plt} B**).
These shared 41 genomes held an additional 10.3% of variable importance over the shared hashes because some genomes contain additional hashes not shared across all models.
Using sourmash lca classify to assign GTDB taxonomy, we find 38 species represented among the 41 genomes (**Figure {@fig:genomes-plt} B**). 

![Some predictive hashes from random forest classifiers anchor to known genomes. **A** 75.1-80.3% of all hashes used to train classifiers anchor to known genomes in RefSeq, GenBank, or human microbiome metagenome-assembled genome databases. A further 4.2-5.6% of hashes anchor to metapangenomes of a subset of these genomes. **B** The 3,859 hashes shared between at least five classifiers anchor to 41 genomes. Genomes account for different amounts of variable importance in each model. Genomes are labelled by 38 GTDB taxonomy assignments. **C** Jaccard similarity between 41 genomes. The highest similarity between genomes is 0.37 and is shared by genomes of the same species, while most genomes have no similarity. This indicates that each genome represents distinct nucleotide sequence.](./figures/fig4.pdf){#fig:genomes-plt}


### Marker genes dominate signatures of IB

While k-mer based signatures allow us to use all sequencing data in a metagenome and quickly compare against all known genomes, hashes lack sequence context and do not represent function. 
Given this, we next sought to uncover the functional potential of the shared hashes. 
To annotate each hash, we reasoned that a hash with predictive importance would be nearby the genomic feature driving the predictive signal in both genomic DNA and the DNA assembly graph. 
Therefore, we performed assembly graph queries on each metagenome using each gene in the 51 shared genomes. 
We then identified all gene query neighborhoods in which each hash occurred, and transferred those annotations to the hash.

We noticed that many hashes were annotated as marker genes: 440 hashes accounting for 7.5% of variable importance across all models annotated to bacterial single copy marker genes [@parks2015checkm; @na2018ubcg], as well as 16S and 23S ribosomal RNA. 
Given the substantial fraction of variable importance attributable to these genetic elements, we were curious how well models built from marker genes alone would perform in IBD subtype classification. 
However, we wanted to only look at marker gene abundances from the 41 shared genomes. 
We first performed cDBG queries using the 41 genomes to retreive all reads in the assembly graph neighborhoods of those genomes. 
We then built random forest classifiers using the same approach as with hashes, but using abundance of 14 ribosomal marker genes and 16S rRNA [@woodcroft2018singlem].
Classification accuracy across all models was similar to the k-mer based model (**Table {@tbl:acc}**), however performs marginally better at CD classification and marginally worse at UC classification (**Figure SUPPLEMENTAL CONFUSION MATRICES**) 
Both model types ranked a 16S rRNA sequence from the genus *Acetatifactor* as having the highest variable importance across studies, demonstrating that while based on different data features, both model types extract similar information.  

|**Validation Study** | **Hash model** | **Ribosomal model** | **Hash model of ribosomal reads** | 
|-----------------|------------|-----------------|------------|
|SRP057027        |    75.9    |      86.4       |   71.7   |
|PRJNA237362      |    71.4    |       75        |   64.3   |
|PRJEB2054        |    69.4    |      19.1       |   15.5   |
|PRJNA385949      |    52.9    |      52.9       |   41.2   |
|PRJNA400072      |    50.9    |      48.1       |   47.4   |
|iHMP             |    49.1    |      44.2       |   46.5   |
Table: Accuracy of random forest classifiers built with different underlying representations of IBD metagenomes when applied to each validation set. {#tbl:acc}

The marker gene model performed similarly as the k-mer model for all studies except PRJEB2054, however this study was sequenced with 36 base pair reads.
It performed poorly due to decreased prediction of marker genes from reads. 
While the hash model performs well even with very short reads, this has limited technological application given ever-increasing sequencing read lengths.

While hash and ribosomal models performed similarly, we were curious whether they captured the same information or whether they captured overlapping but distinct characteristics about IBD subtype. 
Therefore, we used reads used to build the marker gene model to build a hash model from marker genes alone. 
These hash models of ribosomal reads performed worse than both the full hash model and ribosomal model (**Table {@tbl:acc}**}, indicating that each model captures separate, additonal information important for IBD classification.
The ribosomal model likely out-performed the hash model of ribosomal reads given both subsampling and shorter sequence lengths used by k-mers, both of which obscure taxonomic relationships between short sequences.

To test whether marker gene model accuracies were uniquely driven by marker gene abundance in the 41 shared genomes, we next built a model from all marker genes in the metagenome. 
XXXXX.

### Decreased diversity punctuated by strain enrichment explains IBD

To better understand microbial signatures of IBD captured by the hash model, we performed differential abundance analysis on the shared hashes [@martin2020modeling].
Many more hashes were differentially abundant in CD than UC (1815 hashes versus 166 hashes, respectively), reflecting the heterogeneity of the UC gut microbiome [@franzosa2019].
Differentially abundant marker genes segregated to the fraction of hashes that were decreased in abundance in CD. 
This, combined with results from results from marker gene models, indicates taht decrease in diversity captured by marker genes is a hallmark of the CD gut microbiome. 
The majority of non-marker genes with decreased abundance in CD appear to be driven by general decrease in diversity of the gut microbiome, not systematic loss of specific functional potential.

Two strains of *Faecalibacterium prausnitzii* have the largest number of hashes with decreased abundance compared to nonIBD. 
*F. prausnitzii* is an obligate anaerobe and a key butyrate producer in the gut, thereby playing a crucial role in reducing health inflammation [@lopez2017faecalibacterium].
*F. prausnitzii* is extremely sensitive to oxygen, though may be able to withstand exposure for up to 24 hours depending on the availability of metabolites for extracellular electron transfer [@lopez2017faecalibacterium].
Hashes annotated to *Acetatifactor* (GTDB species *s__Acetatifactor sp900066365*) held the largest variable importance of hashes with decreased abundance compared to nonIBD.   
*Acetatifactor* is a bile-acid producing bacteria associated with a health gut, but limited evidence has associated it with decreased abundance in IBD [@pathak2018intestine].
In UC, *Gemmiger formicilis* held both the largest variable importance and number of hashes with decreased abundance compared to non-IBD.
*G. formicilis* is a strictly anaerobic bacteria that produces both formic acid butyric acid [@gossling1975gemmiger].

While most shared genomes decrease in abundance, a substantial portion of hashes from four genomes in CD and two in UC are more abundant in disease. 
These four belong to *Faecalicatena gnavus* (referred to as *[Ruminococcus] gnavus* in NCBI taxonomy and IBD literature) and  *Clostridium bolteae* in CD, and two genomes in the *Flavonifractor* in CD and UC. 
KEGG enrichment analysis on more and less abundant hashes deomonstrated enrichment of pathways dominated by marker genes (e.g. Ribosomes, tRNA biosynthesis) in the less abundant hashes from CD in all four genomes, indicating a general decrease in abundance for these organisms. 
However, metabolic pathways such as starch and sucrose metabolism, propanoate metabolism, and peptidoglycan synthesis (MAKE FIGURE) are enriched in the more abundant hashes in CD.
This shift is indicative of strain-specific enrichment in CD.  

ADD:

+ kmer accumulation curves
+ ...


### Other diff abund bio results

#### c bolt
Given these associations, we performed differential abundance analysis on the *C. bolteae* pangenome between CD and nonIBD.
We compared our results against study of virulence-causing gene in *C. bolteae* [@lozupone2012identifying], and find that 24 of 41 previously identified orthologs are significantly induced in CD. 
Seven of these orthologs are associated with response to oxidative stress. 
(OXIDATIVE STRESS IBD BIO TIE IN). 

We then performed gene enrichment analysis on the differentially abundant genes with KEGG ortholog annotations in *C. bolteae*. 
While many KEGG pathways are significant, flagellar assembly had the second lowest p value (17 genes).
Bacterial flagellin is a dominant antigen in Crohn's disease but not ulcerative colitis [@lodes2004bacterial; @duck2007isolation]. 
#### f gnavus
We performed differential abundance analysis between CD and nonIBD as well as UC and non IBD to understand whether the metapangenome varied between disease states. 
While 5,984 genes were differentially abundant in CD, only 197 were less abundant in UC.
This suggests that *F. gnavus* is different from nonIBD in CD alone. 

We next investigated whether the gene cluster thought to be involved in biosynthesis of the inflammatory polysaccharide was significantly induced in CD. 
We identified 19 of 23 ORFs in the *F. gnavus* pangenome that matched the putative genes in the cluster, all of which were more abundant in CD. 
Further, two subsets, one containing 5 ORFs and one contain 7, were co-located on two contiguous sequences, indicating these genes do form a cluster. 
We then investigated whether this gene cluster was present in non-IBD samples, and found an average of more than 100 reads that mapped per gene in the cluster in 10 of 213 nonIBD metagenomes. 
This indicates that while more abundant in CD, it is also identifiable within healthy human gut microbiomes.    
     
We also genes involved in oxidative stress resistance that are more abundant in CD. 
This includes one super oxide dismutase and five NADH oxidases.

While this evidence supports the idea that *F. gnavus* is harmful in CD, we see some genes that are more abundant in CD that are beneficial for gut health. 
For example, we find 10 a-L-fucosidases. Tryptophan metabolism. ?

### Operons in differentially abundant genes (tmp title)

Given that all genes detected from the *F. gnavus* inflammatory polysaccharide biosynthetic gene cluster were significantly induced in CD, and that subsets of these sequences were colocated on single contiguous sequences, we reasoned that other biologically meaningful genes were likely to occur in clusters. 
Using results from differential abundance analysis, we searched for gene clusters of five or more genes. 
We selected five as a signal:noise compromise, as five was the smallest consecutive stretch detected in the *F. gnavus* cluster. 

We find no evidence of gene clusters that are more abundant in UC. 
Conversely, we find many gene clusters in XX pangenomes that are more abundant in CD. 
XXX



### Predictive hashes not in the metapangenomes XXX

+ *9.1% of hashes*
+ *sgc query by hash*
+ *Assemble, deepvirfinder, mifaser, compare to viral db, etc.*

## Discussion

We present XXX.

In this investigation, we find that gut microbiomes from both UC and CD suffer from stochastic loss of diversity. 

While *C. bolteae* and *R. gnavus* emerge as bad actors in the pathophysiology of CD, no similar signal is detected for UC.
This suggests that while both diseases are associated with lower diversity, CD is uniquely exacerbated by microbes that become more abundant during disease. 

## Methods

All code associated with our analyses is available at www.github.com/dib-lab/2020-ibd/

### IBD metagenome data acquisition and processing

We searched the NCBI Sequence Read Archive and BioProject databases for shotgun metagenome studies that sequenced fecal samples from humans with Crohn's disease, ulcerative colitis, and healthy controls. 
We included studies sequenced on Illumina platforms with paired-end chemistries and with sample libraries that contained greater than one million reads. 
For time series intervention cohorts, we selected the first time point to ensure all metagenomes came from treatment-naive subjects. 

We downloaded metagenomic fastq files from the European Nucleotide Archive using the "fastq_ftp" link and concatenated fastq files annotated as the same library into single files. 
We also downloaded iHMP samples from idbmdb.org.
We used Trimmomatic (version 0.39) to adapter trim reads using all default Trimmomatic paired-end adapter sequences (`ILLUMINACLIP:{inputs/adapters.fa}:2:0:15`) and lightly quality-trimmed the reads (`MINLEN:31 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2`) [@bolger2014]. 
We then removed human DNA using BBMap and a masked version of hg19 [@bushnell2014]. 
Next, we trimmed low-abundance k-mers from sequences with high coverage using khmer's `trim-low-abund.py` [@crusoe2015].  

Using these trimmed reads, we generated scaled MinHash signatures for each library using sourmash (k-size 31, scaled 2000, abundance tracking on) [@brown2016]. 
At a scaled value of 2000, an average of one k-mer will be detected in each 2000 base pair window, and 99.8% of 10,000 base pair windows will have at least one k-mer representative. 
We selected a k-mer size of 31 because of its species-level specificity [@koslicki2016].
A signature is composed of hashes, where each hash represents a k-mer contained in the original sequence.
We retained all hashes that were present in multiple samples, and refer to these as filtered signatures. 

### Principle Coordinates Analysis

We used jaccard distance and cosine distance implemented in `sourmash compare` to pairwise compare filtered signatures. 
We then used the `dist()` function in base R to compute distance matrices. 
We used the `cmdscale()` function to perform principle coordinate analysis [@gower1966]. 
We used ggplot2 and ggMarginal to visualize the principle coordinate analysis [@wickham2019]. 
To test for sources of variation in these distance matrices, we performed PERMANOVA using the `adonis` function in the R vegan package [@oksanen2010].
The PERMANOVA was modeled as `~ diagnosis + study accession + library size + number of hashes`.

### Random forest classifiers

We built random forests classifier to predict CD, UC, and non-IBD status using filtered signatures (hash models), marker genes in the shared 41 genomes (marker gene models), signatures from reads that were detected as marker genes (hash models of marker genes), and marker genes in the full metagenome (full marker gene models).
 
For models from signatures, we transformed sourmash signatures into a hash abundance table where each metagenome was a sample, each hash was a feature, and abundances were recorded for each hash for each sample. 
We normalized abundances by dividing by the total number of hashes in each filtered signature. 
We then used a leave-one-study-out validation approach where we trained six models, each of which was trained on five studies and validated on the sixth. 
To build each model, we first performed vita variable selection on the training set as implemented in the Pomona and ranger packages [@degenhardt2017; @wright2015]. 
Vita variable selection reduces the number of variables (e.g. hashes) to a smaller set of predictive variables through selection of variables with high cross-validated permutation variable importace [@janitza2018].
Using this smaller set of hashes, we then built an optimized random forest model using tuneRanger [@probst2019]. 
We evaluated each validation set using the optimal model, and extracted variable importance measures for each hash for subsequent analysis. 
To make variable importance measures comparable across models, we normalized importance to 1 by dividing variable importance by the total number of hashes in a model and the total number of models.  

For the marker gene models, we generated marker gene abundances for 14 ribosomal marker genes and 16S rRNA using singleM [@woodcroft2018singlem].
We then followed the same model building procedure as the hash models. 
 
### Anchoring predictive hashes to genomes

We used sourmash `gather` with parameters `k 31` and `--scaled 2000` to anchor predictive hashes to known genomes [@brown2016]. 
Sourmash `gather` searches a database of known k-mers for matches with a query [@pierce2019].
We used the sourmash GenBank database (2018.03.29, https://osf.io/snphy/), and built three additional databases from medium- and high-quality metagenome-assembled genomes from three human microbiome metagenome reanalysis efforts (https://osf.io/hza89/) [@pasolli2019; @nayfach2019; @almeida2019].
In total, approximately 420,000 microbial genomes and metagenome-assembled genomes were represented by these four databases. 
We used the sourmash `lca` commands against the GTDB taxonomy database to taxonomically classify the genomes that contained predictive hashes.
To calculate the cumulative variable importance attributable to a single genome, we used an iterative winner-takes-all approach.
The genome with the largest fraction of predictive k-mers won the variable importance for all hashes contained within its genome.
These hashes were then removed, and we repeated the process for the genome with the next largest fraction of predictive k-mers. 

To identify hashes that were predictive in at least five of six models, we took the union of predictive hashes from all combinations of five models, as well as from the union of all six models.
We refer to these hashes as shared predictive hashes.
We anchored variable importance of these shared predictive hashes to known genomes using sourmash `gather` as above. 

### Compact de Bruijn graph queries for predictive genes and genomes

To annotate hashes with functional potential, we first extracted open reading frames (ORFs) from the shared 41 genomes using prokka, and annotated ORFs with EggNog [@seemann2014prokka; @huerta2019eggnog].
When then used spacegraphcats `multifasta_query` to create a hash:gene map.
Spacegraphcats retreives k-mers in the compact de Bruijn graph neighborhood of a query gene, and hashing these k-mers via sourmash generates a hash:gene map [@brown2020exploring; @brown2016]. 
Because genomes with shared 31-mers may annotate the same hash, we allowed hashes to be annotated multiple times. 
This was particularly appropriate for hashes from highly conserved regions, e.g. 16S ribosomal RNA.

We used spacegraphcats `search` to retreive k-mers in the compact de Bruijn graph neighborhood of the shared genomes [@brown2020exploring]. 
We then used spacegraphcats `extract_reads` to retreive the reads and `extract_contigs` to retreive unitigs in the compact de Bruijn graph that contained those k-mers, respectively.
These reads were used to generate marker gene abundances for the 41 shared genomes for the marker gene random forest models.
 
### Differential hash abundance analysis

To determine whether shared hashes were differentially abundant from nonIBD in UC or CD, we used corncob [@martin2020modeling].
We used all hash abundances from sourmash signatures to determine hash library size, and then compared hash abundances between disease groups using the likelihood ratio test with the formula `study_accession + diagnosis` and the null formula `study_accession` [@martin2020modeling]. 
We considered genes with p values < .05 after bonferonni correction as statistically significant. 
We performed enrichment analysis using the R package clusterProfiler [@yu2012clusterprofiler]. 

## References
