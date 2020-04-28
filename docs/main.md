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

Here we perform a meta-analysis of six studies of IBD gut metagenome cohorts comprising 260 CD, 132 UC and 213 healthy controls [@lloyd2019; @lewis2015; @hall2017; @franzosa2019; @gevers2014; @qin2010]. 
First, we re-analyzed each study using a consistent k-mer-based, reference-free approach. 
We demonstrate that study and patient predict more variation than disease status. 
Next, we used random forests to accurately predict disease status and to determine the k-mers that are predictive of UC and CD. 
Then, we use compact de Bruijn graph queries to reassociate k-mers with sequence context and perform taxonomic and functional characterization of these sequence neighborhoods. 
Our analysis pipeline is lightweight and relies on well-documented and maintained software, making it extensible to other large cohorts of metagenomic sequencing data.

## Results

### Many cohorts of IBD stool metagenome samples

Sample collection, DNA extraction, library preparation, and sequencing all introduce technical artifacts and biases into metagenomics sequencing data that are difficult to separate from biological signal. 
Metaanlysis improves detection of biological signal by obscuring signals of technical artifacts that vary between studies. 
As such, we used six IBD cohorts with stool metagenomes in this analysis {@tbl:cohorts}. 
We selected cohorts with Illumina shotgun stool metagenomes for which CD, UC, and nonIBD status was recorded for all samples.
Any sample with fewer than one million reads was eliminated. 
All studies with time-series data were used as validation sets. 
In total, we analyzed 605 samples.   


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

### Reference-free approach for meta-analysis of IBD metagenomes.

Given that both reference-based and *de novo* methods suffer from substantial and biased loss of information in the analysis of metagenomes, we sought a reference- and assembly- free pipeline to fully characterize each sample. 
We used scaled MinHash sketching as implemented in sourmash to produce a compressed representation of the k-mers contained in each sample [@pierce2019]. 
K-mers are nucleotide sequences of length *k*. 
K-mers are superior to alignment and assembly in metagenome analysis because: 
1) k-mers enable exact matching, which is fast and requires little computational resources; 
2) k-mers do not need to be present in reference databases to be included in analysis; and 
3) k-mers capture information from reads even when there is low coverage or high strain variation, both of which preclude assembly. 
We used a k-mer size of 31 in our analysis as high similarity (e.g. many k-mers in common) at this size approximates species-level taxonomic similarity [@koslicki2016]. 
Further, a 4.7 Mbp genome such as *Escherichia coli* K-12 (substrain MGI655) contains approximately 4.6 million k-mers, but only a fraction of these are needed to recapitulate similarity measurements using other metrics [@pierce2019]. 
Therefore, we retained 1/2000th of the distinct, random k-mers from each sample or our analysis. 
At a scaled value of 2000, an average of one k-mer will be detected in each 2000 base pair window, and 99.8% of 10,000 base pair windows will have at least one k-mer representative. 
Because the same sets of k-mers are retained, this approach creates a consistent set of marker sequences across samples. 
This enables comparisons between samples.

Using this method, adapter sequences, human DNA, and sequencing errors can falsely inflate or deflate similarity measurements. 
As such, we also used a consistent preprocessing pipeline to adapter trim, remove human DNA, and k-mer trim (e.g. remove erroneous k-mers). 
Because k-mer trimming retains some erroneous k-mers, we further filtered signatures to remove hashes that occured once across all samples. 
This removed additional k-mers that were most likely to be errors while keeping k-mers that were real but low abundant in some samples. 
There were 57,479,783 distinct hashes across all samples, of which 9,334,204 remained after filtering. 


[//]: # The most common hash appeared in 944 of 954 non iHMP samples.

### Study accounts for more variation that disease across cohorts

In this study, we aimed to identify microbial signatures associated with IBD. 
However, given that biological and technical artifacts can differ greatly between metagenome studies [@wirbel2019], we first quantified these sources of variation. 
Only variables patient, study, diagnosis, and library size were recorded for all samples. 
Using filtered signatures, we calculated pairwise distance matrices using jaccard distance and cosine distance, where jaccard distance captured sample richness and cosine distance captured samle diversity. 
We performed principle coordinate analysis and PERMANOVA using these distance matrices (**Figure {@fig:comp-plts}**). 
Patient accounts for the most variation of the four factors tested (jaccard: .622, p < .001; cosine: .647, p < .001). 
This reinforces that inter-individual variation is higher than intra-individual variation (CITATIONS). 
Study accounts for the next highest amount of variation (jaccard: .122, p < .001; cosine: .118, p < .001). 
This illustrates that technical artifacts introduce biases with strong signals that can be ameliorated by metanalysis. 
Diagnosis accounts for a small but significant amount of variation (jaccard: .035, p < .001; cosine .020, p < .001), a pattern we would expect given the heterogeneity of the IBD disease spectrum and echoes patterns observed for colorectral cancer [@wirbel2019]. 
We conclude that both sequence richness and diversity vary by patient, study, and disease status.

![Principle component analysis of metagenomes from IBD cohorts. 
Individual and study account for the most variation in for both A) Jaccard distance and B) PCoA of cosine distance.](../figures/minhash-comp-plts-all.pdf){#fig:comp-plts}
  
### Metagenomic markers differentiate IBD status across cohorts

To evaluate whether filtered signatures predict IBD status across studies, we built a random forests classifier. 
We used 70% of metagenome samples from three non-time series cohorts as a training set, the remaining 30% as a test set, and three time series cohorts as validation sets. 
Given the high-dimensional structure of this dataset (e.g. many more k-mers than samples), we first used the vita method to select predictive hashes in the training set [@degenhardt2017]. 
This reduced the number of hashes used in the classifier from 9,334,204 to 14,481. 
Vita variable selection iteratively removes the bottom 20% of predictive k-mers in successive rounds of random forests, and stops when the out of bag error begins to rise. 
In omics data, this results in the elimination of correlated features that carry redundant information [@he2010]. 

Using the predictive hashes, we then built a random forest classifier to predict UC, CD, or non-IBD status across samples. 
100% of samples were classified accurately in the training set  while 92% of samples were classified accurately in the test set. 

We validated our model on three time series cohorts. 
Given that patient explains the most variation in samples, we excluded these cohorts from the training dataset to reduce bias from intra-individual testing and training. 
Accuracy ranged from 65.8%-80.4% in our validation cohorts (**Table ?**). 
Interestingly, out model was robust to read preprocessing techniques; in cohorts SRP057027 and PRJNA385949, we used the same read preprocessing pipeline as we used for our training and testing sets. 
Accuracy was 80.4% and 65.4% for these cohorts, respectively. 
Then, for the iHMP cohort, we used preprocessed reads as provided at ibdmdb.org, which used a different read quality control pipeline.
Accuracy was 77.6% on this cohort.

[//]: # SRP057027; 337 samples, 113 patients) @lewis2015 all CD

[//]: # PRJNA385949; 155 samples, 24 patients) @hall2017

[//]: # iHMP IBD cohort (1338 samples, 106 patients) 

The most common error our model made was misclassification of CD as nonIBD. 
We next sought to understand the source of this error. 
CD presents mouth-anus and can localize to the upper gastrointestinal tract, the illeum, the colon, or any combination of those sites. 
We hypothesized that CD misclassification as nonIBD would occur more frequently in non-colonic CD patients. 
The iHMP IBD cohort has extensive metadata, including baseline Montreal classification for disease location for most CD patients. 
Using this localization information, we found that non-colonic CD patients were significantly more likely to be misclassified as nonIBD than colonic patients (chi-square test, p < .001). 


[//]: # compare to other random forests for IBD from literature

### Some predictors of IBD associate with known genomes

We next evaluated the identity and function of the predictive k-mers in our random forests classifier. 
We first compared the k-mers against sequences in reference databases. 
We used sourmash gather to estimate the proportion of sequenced genomes contained in our predictive k-mers [@pierce2019]. 
We compared our predictive hashes against all microbial genomes in GenBank, as well as metagenome-assembled genomes from three recent reassembly efforts from human microbiome metagenomes [@pasolli2019; @nayfach2019; @almeida2019]. 
68.7% of hashes matched 129 genomes from these databases. 
The largest number of hashes matched *Clostriales bacterium* (2.9% hashes, 18.2% of genome), *Ruminococcus gnavus* (2.8% hashes, 21.2% of genome), and *Clostridium clostridioforme* (2.7% of hashes, 13.5% of genome) (**Figure {@fig:gather}**). 
These genomes have been associated with IBD in past studies (CITATIONS).


![Taxonomic composition of predictive k-mers present in reference databases. 
We used sourmash gather to identify genome matches in GenBank and human microbiome metagenome-assembled genomes. 
Then, we used GTDBtk to assign taxonomy to each match.](../figures/gather_order_viz.pdf){#fig:gather}

### Unknown but predictive k-mers XXX

Given that 31.3% of predictive k-mers from our random forest classifier did not match known sequences in databases, we next sought to characterize these k-mers. 
We reasoned that many unknown but predictive k-mers likely originate from closely related strain variants of identified genomes, thus we first sought to recover these associations. 
We performed compact de Bruijn graph queries with the 129 genomes that contained predictive k-mers (sgc CITATION?) to extract neighboring sequences to these genomes. 
85.4% of unknown k-mers were in the neighborhood of these 129 genomes (68.7% in query genomes, 16.7% in query neighborhoods).
This suggests that 16.7% of predictive k-mer originate from strain-variable or accessory genome components. 
These components are not recoverable by reference-based or *de novo* approaches, but are important for disease classification.

An additional 14.6% of k-mers were neither in known genomes or in the neighborhood of known genomes. ... 

## Discussion

+ This pipeline generates an inclusive summary of all sequences contained in the metagenome while reducing the computational footprint of the data and its analysis.
+ (CD misclassification results) This indicates that the colon microbiome of patients with non-colonic CD presents more like a microbiome of a nonIBD patient, suggesting that microbial signatures of IBD are localized to the site of disease. 
+ Our random forest results indicate that there are markers of IBD that are conserved across studies. 
 
## Methods

All code associated with our analyses is available at www.github.com/dib-lab/2020-ibd/

### IBD metagenome data acquisition and processing

We searched the NCBI Sequence Read Archive and BioProject databases for shotgun metagenome studies that sequenced fecal samples from humans with Crohn's disease, ulcerative colitis, and healthy controls. 
We included studies sequenced on Illumina platforms with paired-end chemistries and with sample libraries that contained greater than one million reads. 

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
