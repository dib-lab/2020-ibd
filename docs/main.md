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
geometry: margin=1in
bibliography: bibliography.bib
header-includes:
    - \usepackage{setspace}
    - \doublespacing
    - \usepackage{lineno}
    - \linenumbers
fontfamily: helvet
---

## Introduction

Metagenomic sequencing captures the functional potential of microbial communities through DNA sequencing of genes and organisms. 
Metagenomics has been used to profile many human microbial communities, including those that change in or contribute toward disease. 
In particular, human gut microbiomes have been extensively characterized for their potential role in diseases such as obesity [@greenblum2012], type II diabetes [@qin2012], colorectal cancer [@wirbel2019], and inflammatory bowel disease [@lloyd2019; @morgan2012; @hall2017; @franzosa2019]. 
Inflammatory bowel disease (IBD) refers to a spectrum of diseases characterized by chronic inflammation of the intestines and is likely caused by host-mediated inflammatory responses at least in part elicited by microorganisms [@kostic2014]. 
However, no causative or consistent microbial signature has been associated with IBD to date. 

[//]: # need a transition sentence here

Statements about biology...

[//]: # In many studies, IBD is associated with lower biodiversity [@lewis2015; @moustafa2018; @qin2010], however one study failed to replicate this finding [@lloyd2019]. Some studies have noted lower Firmicutes and higher Gammaproteobacteria in IBD, however one study noted that disease state is not independently associated with microbiome shifts after correcting for other covariates [@morgan2012]. Another study found that a novel *Ruminococcus gnavus* clade is enriched in IBD and that this enrichment would be masked in 16s rRNA amplicon studies given commensurate shifts in other closely-related species [@hall2017]. Additionally, time-course studies of IBD have shown that inter-individual differences contribute more to microbiome variation than disease status [@lloyd2019]. 

[//]: # Both reactive oxygen species and reactive nitrogen species are likely important in IBD [@keshavarzian2003; @winter2014]. 

Although there is no consistent taxonomic or functional trend in the gut microbiome associated with IBD status, metagenomic studies conducted unto this point have left substantial portions of reads unanalyzed. 
Reference-based pipelines commonly used to analyze metagenomic data from IBD cohorts such as HUMANn2 characterize on average 31%-60% of reads from the human gut microbiome [@franzosa2014; @lloyd2019]. 
Reads fail to map when there is no closely related organism or sequence in the reference database. 
Reads that do not map to references are typically ignored in downstream analysis. 
To combat these issues, reference-free approaches like *de novo* assembly and binning are used to generate metagenome-assembled genome bins (MAGs). 
MAGs represent species-level composites of genomes of organisms in a sample, and thus often more closely recapitulate the organisms in the sample. 
However, *de novo* approaches fail when there is low-coverage of or high strain variation in gut microbes, or with sequencing error [@olson2017]. 
Even when performed on a massive scale,an average of 12.5% of reads fail to map to all *de novo* assembled organisms from human microbiomes [@pasolli2019], meaning some sequences are not assembled or binned.
As with reference-based approaches, these reads are typically left unanalyzed.

Further, within-study confounders from biological and technical artifact may lead to false associations. 
This can be ameliorated by meta-analysis.  

Here we perform a meta-analysis of six studies of IBD gut metagenome cohorts comprising 260 CD, 132 UC and 213 healthy controls (see **Table {@tbl:cohorts}**) [@lloyd2019; @lewis2015; @hall2017; @franzosa2019; @gevers2014; @qin2010]. 
First, we re-analyzed each study using a consistent k-mer-based, reference-free approach. 
We demonstrate that study and patient predict more variation than disease status. 
Next, we used random forests to accurately predict disease status and to determine the k-mers that are predictive of UC and CD. 
Then, we use compact de Bruijn graph queries to reassociate k-mers with sequence context and perform taxonomic and functional characterization of these sequence neighborhoods. 
Our analysis pipeline is lightweight and relies on well-documented and maintained software, making it extensible to other large cohorts of metagenomic sequencing data.

## Results

| Cohort      | Cohort names             | Country          | Total | CD  | UC  | nonIBD | Reference    |
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

Given that both reference-based and *de novo* methods suffer from substantial and biased loss of information in the analysis of metagenomes, we sought a reference- and assembly-free pipeline to fully characterize each sample.
K-mers, words of length *k* in nucleotide sequences, have previously been exploited for annotation-free characterization of sequencing data (CITATION).
K-mers are superior to alignment and assembly in metagenome analysis because: 
1) k-mers enable exact matching, which is fast and requires little computational resources; 
2) k-mers do not need to be present in reference databases to be included in analysis; and 
3) k-mers capture information from reads even when there is low coverage or high strain variation, both of which preclude assembly. 
However, k-mers are complex sets given that there are approximately n^k k-mers in a nucleotide sequence (e.g. 4.7 Mbp genome such as *Escherichia coli* K-12 (substrain MGI655) contains approximately 4.6 million k-mers).  
However, only a fraction of these are needed to recapitulate similarity measurements using other metrics [@pierce2019]. 
Thus we used scaled MinHash sketching as implemented in sourmash to produce a compressed representation of the k-mers contained in each sample [@pierce2019]. 
At a scaled value of 2000, an average of one k-mer will be detected in each 2000 base pair window, and 99.8% of 10,000 base pair windows will have at least one k-mer representative. 
We refer to the subsampled representative set of k-mers as a *signature*, and to each subsampled k-mer in a signature as a *hash*. 
Importantly, this approach creates a consistent set of hashes across samples by retaining the same hashes when the same k-mers are observed. 
This enables comparisons between metagenomes.

Using this method, adapter sequences, human DNA, and sequencing errors can falsely inflate or deflate similarity measurements. 
As such, we also used a consistent preprocessing pipeline to adapter trim, remove human DNA, and k-mer trim (e.g. remove erroneous k-mers). 
Because k-mer trimming retains some erroneous k-mers, we further filtered signatures to retain hashes that were present in multiple signatures. 
This removed hashes that were likely to be errors while keeping hashes that were real but of low abundance in some signatures. 
There were 46,267,678 distinct hashes across all samples, of which 7,376,151 remained after filtering. 

### K-mers capture variation due to disease subtype

In this study, we aimed to identify microbial signatures associated with IBD. 
However, given that biological and technical artifacts can differ greatly between metagenome studies [@wirbel2019], we first quantified these sources of variation. 
We calculated pairwise distance matrices using jaccard distance and cosine distance between filtered signatures, where jaccard distance captured sample richness and cosine distance captured sample diversity. 
We performed principle coordinate analysis and PERMANOVA with these distance matrices (**Figure {@fig:comp-plts}**), using the variables study accession, diagnosis, library size, and number of hashes in a filtered signature (**Table {@tbl:permanova}**). 
Number of hashes in a filtered signature accounts for the highest variation, possibly reflecting reduced diversity in stool metagenomes of CD and UC patients (CITATIONS). 
Study accounts for the second highest variation, emphasizing that technical artifacts can introduce biases with strong signals.
Diagnosis accounts similar amount of variation as study, demonstrating that there is a small but detectable signal of IBD subtype in stool metagenomes.   
We conclude that both sequence richness and diversity vary by number of hashes, study, and diagnosis.

![Principle coordinate analysis of metagenomes from IBD cohorts.](../figures/minhash-comp-plts-all.pdf){#fig:comp-plts}
  
### Hashes are weakly predictive of IBD subtype

To evaluate whether the variation captured by diagnosis is predictive of IBD disease subtype, we built a random forests classifier to predict CD, UC, or non-IBD.
We used a leave-one-study-out cross-validation approach where we built and optimized a classifier using five cohorts and validated on the sixth.  
Given the high-dimensional structure of this dataset (e.g. many more hashes than samples), we first used the vita method to select predictive hashes in the training set [@janitza2018, @degenhardt2017]. 
Vita variable selection is based on permuation of variable importance, where p-values for variable importance are calculated against a null distribution that is built from variables that are estimated as non-important [@janitza2018].
This approach retains important variables that are correlated [@janitza2018; @seifert2019], which is desirable in omics-settings where correlated features are often involved in a coordinated biological response, e.g. part of the same operon, genome, or pathway (CITATIONS). 
Variable selection reduced the number of hashes used in each model to XX-XX (**TABLE?**). 
Using this reduced set of hashes, we then optimized each random forests classifier on the training set, producing six optimized models.
We validated each model on the left-out study.
The accuracy on the validation studies ranged from 49.1%-75.9% (**FIGURE ACCURACY BARPLOT AND CONFUSION MATRICES**), outperforming previously published models built on metagenomic data alone (CITATIONS).

We next sought to understand whether there was a consistent biological signal captured among classifiers by evaluating the fraction of shared hashes between models.
We intersected each set of hashes used to build each optimized classifier (**FIGURE UPSET PLOT WITH VAR IMP**).
Nine hundred thrity two hashes were shared between all classifiers, while 3,859 hashes were shared between at least five studies.
The presence of shared hashes between classifiers indicates that there is a weak but consistent biological signal for IBD subtype between cohorts.      

Shared hashes accounted for XX% of all hashes used to build the optimized classifiers.
If shared hashes are predictive of IBD subtype, we would expect that these hashes would account for an outsized proportion of variable importance in the optimized classifiers.
To calculate the relative variable importance contributed by each hash, we first normalized the variable importance values within each classifier by dividing by the total variable importance (e.g. sum to 1 within each classifier).
We then normalized the variable importance of across all classifiers by dividing by the total number of classifiers (e.g. divided by six so the total variable importance of all hashes across all classifiers summed to 1).
40.2% of the total variable importance was held by the 3,859 hashes shared between at least five classifiers, with XX% attributable to the 954 hashes shared between all six classifiers. 
This indicates that shared hashes contribute a large fraction of predictive power for classification of IBD subtype. 

### Some predictive hashes anchor to known genomes

We next evaluated the identity of the predictive hashes in each classifier. 
We first compared the predictive hashes against sequences in reference databases. 
We used sourmash gather to anchor predictive hashes to known genomes [@pierce2019]. 
We compared our predictive hashes against all microbial genomes in GenBank, as well as metagenome-assembled genomes from three recent reassembly efforts from human microbiome metagenomes [@pasolli2019; @nayfach2019; @almeida2019]. 
Between 75.1-80.3% of of hashes anchored to 1,161 genomes (**FIGURE STACKED BARPLOT DB %**). 
However, the 3,859 hashes shared between at least five classifiers anchored to only 41 genomes (**FIGURE GENOME UPSET PLOT or CUMULATIVE VARIABLE IMPORTANCE PLOT FOR GENOMES**).
Futher, these 41 genomes accounted for 50.5% of the total variable importance, a 10.3% increase over the hashes alone.
In contrast to all hashes, only 69.4% of these hashes were identifiable, a decrease of 5.7-10.9%. 
This indicates that hashes that are more likely to be important for IBD subtype classification are less likely to be anchored to genomes in reference databases.

Using sourmash lca classify to assign GTDB taxonomy, we find 38 species represented among the 41 genomes. 
The genome that anchors the most variable importance is **Acetatifactor sp900066565**. 
(Add %phyla/etc? Is it even worth analyzing these that much when everything changes after spacegraphcats?)
However, we observe that while most genomes assign to one species, 19 assign to one or more distantly related genomes. 
When we take the Jaccard index of these 41 genomes, we observe little similarity despite contamination (**FIGURE SOURMASH CONTAM 41 GENOMES**). 
Therefore, we proceeded with analysis with the idea that each of the 41 genomes is a self-contained entity that captures distinct biology.

### Unknown but predictive hashes represent novel pangenomic elements

Given that 30.6% of hashes shared between at least five classifiers did not anchor to genomes in databases, we next sought to characterize these hashes. 
We reasoned that many unknown but predictive hashes likely originate from closely related strain variants of identified genomes and sought to recover these variants. 
We performed compact de Bruijn graph queries into each metagenome sample with the 41 genomes that contained predictive hashes (CITATION: SPACEGRAPHCATS).
This produced pangenome neighborhoods for each of the 41 genomes.
86.1% of unknown hashes shared between at least five classifiers were in the pangenomes of the 41 genomes, a 16.7% increase over the 41 genomes alone. 
This suggests that at least 16.7% of these hashes originate from novel strain-variable or accessory elements in the pangenomes. 
These components are not recoverable by reference-based or *de novo* approaches, but are important for disease classification.
Further, these pangenomes captured an additional 4.2-5.2% of all predictive hashes from each classifier (**FIGURE % BARPLOT**).
The pangenomes also captured 74.5% of all variable importance, a 24% increase over the 41 genomes alone. 
This indicates that pangenomic variation contributes substantial predictive power toward IBD subtype classification.

Pangenomic neighborhood queries disproportionately impact the variable importance anchored by specific genomes (**FIGURE SCATTER PLOT**).
While most genomes maintained a similar proportion of importance with or without pangenome queries, three pangenomes shifted dramatically.
While an *Acetatifactor* species anchored the most importance prior to pangenome construction, the specific species of *Acetatifactor* switched from *sp900066565*, to *sp900066365*. 
This suggests that pangenome queries might give a more complete picture of the strains involved in IBD (DOES THIS SUGGEST SOMETHING DIFFERENT/BETTER?).   

Conversely, *Faecalibacterium prausnitzii_D* increased from anchoring ~2.9% to ~10.5% of the total variable importance, indicating that substantial pangenomic elements were hidden in for this organism in particular. 
Or something. 

### Differential Abundance of Pangenomes

TBD   

## Discussion

+ This pipeline generates an inclusive summary of all sequences contained in the metagenome while reducing the computational footprint of the data and its analysis.
 
## Methods

All code associated with our analyses is available at www.github.com/dib-lab/2020-ibd/

### IBD metagenome data acquisition and processing

We searched the NCBI Sequence Read Archive and BioProject databases for shotgun metagenome studies that sequenced fecal samples from humans with Crohn's disease, ulcerative colitis, and healthy controls. 
We included studies sequenced on Illumina platforms with paired-end chemistries and with sample libraries that contained greater than one million reads. 
For time series intervention cohorts, we selected the first time point to ensure all metagenomes came from treatment-naive subjects. 

We downloaded metagenomic fastq files from the European Nucleotide Archive using the "fastq_ftp" link and concatenated fastq files annotated as the same library into single files. 
We also downloaded iHMP samples from idbudb.org.
We used Trimmomatic (version 0.39) to adapter trim reads using all default Trimmomatic paired-end adapter sequences (`ILLUMINACLIP:{inputs/adapters.fa}:2:0:15`) and lightly quality-trimmed the reads (`MINLEN:31 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2`) [@bolger2014]. 
We then removed human DNA using BBMap and a masked version of hg19 [@bushnell2014]. 
Next, we trimmed low-abundance k-mers from sequences that have high coverage using khmer's `trim-low-abund.py` [@crusoe2015].  

Using these trimmed reads, we generated scaled MinHash signatures for each library using sourmash (k-size 31, scaled 2000, abundance tracking on) [@brown2016]. 
We selected a k-mer size of 31 because of its species-level specificity [@koslicki2016].
A signature is composed of hashes, where each hash represents a k-mer contained in the original sequence; hashing k-mers to integers reduces storage space and improves computational run times when performing comparisons.

Although we adapter, quality, and k-mer trimmed our reads, some erroneous k-mers were likely included in the MinHash sequences, especially for low-coverage sequences. 
Therefore, we retained all hashes that were present in multiple samples.
We refer to these as filtered signatures. 

### Principle Coordinates Analysis

We used jaccard distance and cosine distance implemented in `sourmash compare` to pairwise compare filtered signatures. 
We then used the `dist()` function in base R to compute distance matrices. 
We used the `cmdscale()` function to perform principle coordinate analysis [@gower1966]. 
We used ggplot2 and ggMarginal to visualize the principle coordinate analysis [@wickham2019]. 
To test for sources of variation in these distance matrices, we performed PERMANOVA using the `adonis` function in the R vegan package [@oksanen2010].
The PERMANOVA was modeled as `~ diagnosis + study accession + library size + number of hashes`.

### Random forest classification

We built a random forest classifier to predict CD, UC, and non-IBD status using filtered signatures. 
First, we transformed sourmash signatures into a hash abundance table where each metagenome was a sample, each hash was a feature, and abundances were recorded for each hash for each sample. 
We normalized abundances by dividing by the total number of hashes in each filtered signature. 
We then used a leave-one-study-out validation approach where we trained six models, each of which was trained on five studies and validated on the sixth. 
To build each model, we first performed vita variable selection on the training set as implemented in the Pomona and ranger packages [@degenhardt2017; @wright2015]. 
Vita variable selection reduces the number of variables (e.g. hashes) a smaller set of predictive variables through selection of variables with high cross-validated permutation variable importace [@janitza2018].
Using this smaller set of hashes, we then built an optimized random forest model using tuneRanger [@probst2018]. 
We evaluated each validation set using the optimal model, and extracted variable importance measures for each hash for subsequent analysis. 

[//]: # To test whether CD misclassifications were more common in non-colonic IBD, we used the iHMP metadata (available at ibdmdb.org) to assess baseline Montreal location as defined previously [@silverberg2005]. 
[//]: # We defined L1, L4, and L1+L4 as non-colonic presentation of IBD, and designated all other classification as colonic. 
[//]: # We then used a chi-squared test using the R `chi.square()` function to test whether we observed more misclassification for non-colonic CD patients than for colonic CD patients. 

### Characterization of predictive k-mers

We used sourmash `gather` with parameters `k 31` and `--scaled 2000` to anchor predictive hashes to known genomes [@brown2016]. 
Sourmash `gather` searches a database of known k-mers for matches with a query [@pierce2019].
We used the sourmash GenBank database (2018.03.29, https://osf.io/snphy/), and built three additional databases from medium- and high-quality metagenome-assembled genomes from three human microbiome metagenome reanalysis efforts (https://osf.io/hza89/) [@pasolli2019; @nayfach2019; @almeida2019].
In total, approximately 420,000 microbial genomes and metagenome-assembled genomes were represented by these four databases. 
We used the sourmash `lca` commands against the GTDB taxonomy database to taxonomically classify the genomes that contained predictive hashes.
To calculate the cumulative variable importance attributable to a single genome, we used an iterative winner-takes-all approach.
The genome with the largest fraction of predictive k-mers won the variable importance for all hashes contained within its genome.
These hashes were then removed, and we repeated the process for the genome with the next largest fraction of predictive k-mers. 

To identify hashes that were predictive in at least five of six models, we took the union of predictive hashes from all combinations of five models, as well as from the union of all six models.
We re-anchored these hashes to known genomes using sourmash `gather` as above. 

### Compact de Bruijn graph queries for predictive genomes

We used spacegraphcats `search` to retreive k-mers in the compact de Bruijn graph neighborhood of the genomes that matched predictive k-mers (CITATION). 
We then used spacegraphcats `extract_reads` to retreive the reads and `extract_contigs` to retreive unitigs in the compact de Bruijn graph that contained those k-mers, respectively.

### Characterization of graph pangenomes

**Pangenome signatures** To evaluate 

**Differential abundance**

### Compact de Bruijn graph queries for unknown predictive k-mers

We used the spacegraphcats `query_by_hashval` with a radius of 5 to retreive the compact de Bruijn graph neighborhood of unknown predictive k-mers.

## References
