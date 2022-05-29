# BiocManager::install("WGCNA") 

library(WGCNA)
library(readr)
library(dplyr)

# input & cleaning ----------------------------------------------------------------
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Read in protein abundance data
counts <- read_tsv("~/Downloads/sig_ccs/tmp_sig_prot_counts_all_for_wgcna.tsv")

# remove aux data cols; transpose so genes are columns
datExpr0 = as.data.frame(t(counts[ , -c(1:2)]))
names(datExpr0) = paste(counts$genome, counts$protein, sep = "_")
rownames(datExpr0) = colnames(counts[-c(1:2)])

# check for how many samples are "missing"

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK # if true, move forward

# if the above statement is false, remove offending genes

# if (!gsg$allOK) {
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0) 
#     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0) 
#     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
# }

#  cluster samples to remove any obvious outliers 

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,12)

pdf(file = "~/Downloads/wgcna_sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#dev.off()

# remove all samples above cut height 15, where the outliers are.
# Plot a line to show the cut
abline(h = 15, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 0 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
# datExpr <- datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#  read in trait data
traitData <- read_tsv("~/github/2020-ibd/inputs/working_metadata.tsv")%>%
  select(library_name, diagnosis) %>%
  distinct()

allTraits <- as.data.frame(traitData)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
samples = rownames(datExpr)
traitRows = match(samples, allTraits$library_name)
datTraits = as.data.frame(allTraits[traitRows, -1])
rownames(datTraits) = allTraits[traitRows, 1]
colnames(datTraits) <- c("diagnosis")
datTraits$diagnosis <- gsub("nonIBD", "0", datTraits$diagnosis)
datTraits$diagnosis <- gsub("UC", "1", datTraits$diagnosis)
datTraits$diagnosis <- gsub("CD", "2", datTraits$diagnosis)
datTraits$diagnosis <- as.numeric(datTraits$diagnosis)
collectGarbage()

# visualize how the clinical traits relate to the sample dendrogram
# re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file = "~/Downloads/wgcna_tmp_dendro_colors.pdf", height = 8, width = 10)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# white means a low value, red a high value, and grey a missing entry.

# save data.
save(datExpr, datTraits, file = "~/Downloads/wgcna_dataInput.RData")


# Network construction and module detection ------------------------------

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "~/Downloads/wgcna_dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# build network

net <- blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "~/Downloads/wgcna_TOM", 
                       verbose = 3)

#  Code chunk 4

# how many modules were identified and what the module sizes are.
# label 0 is for the genes that are outside of all modules
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#  Code chunk 5

# save the module assignment and module eigengene information
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "~/Downloads/wgcna_networkConstruction-auto.RData")


save(net, file = "~/Downloads/wgcna_net.RData")

geneInfo0 = data.frame(gene = names(datExpr),
                       moduleColor = moduleColors)
write.csv(geneInfo0, "~/Downloads/wgcna_protein_modules.csv", quote = F, row.names = F)
# Relating modules to external traits and identifying important genes --------

#  Code chunk 1

# # Load the expression and trait data saved in the first part
# lnames = load(file = "~/Downloads/wgcna_dataInput.RData")
# #The variable lnames contains the names of loaded variables.
# lnames
# # Load network data saved in the second part.
# lnames = load(file = "~/Downloads/wgcna_networkConstruction-auto.RData");
# lnames
# 
# # Define numbers of genes and samples
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# # Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
# moduleTraitCor = cor(MEs, datTraits, use = "p")
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 
# 
# #  Code chunk 3
# 
# sizeGrWindow(20,20)
# # Will display correlations and their p-values
# textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
#                     signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# #pdf(file = "sandbox/wgcna/day1_labeldHeatmap.pdf", height = 8, width = 10)
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(200),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# #dev.off()
# #  Code chunk 4
# 
# # Define variable yan containing the initial_yan column of datTrait
# yan = as.data.frame(datTraits$initial_yan)
# names(yan) = "yan"
# # names (colors) of the modules
# modNames = substring(names(MEs), 3)
# 
# geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# 
# names(geneModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# 
# geneTraitSignificance = as.data.frame(cor(datExpr, yan, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
# 
# names(geneTraitSignificance) = paste("GS.", names(yan), sep="");
# names(GSPvalue) = paste("p.GS.", names(yan), sep="");
# 
# 
# #  Code chunk 5
# 
# module = "green"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for yan",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# 
# #  Code chunk 6
# 
# names(datExpr)
# 
# 
# #  Code chunk 7
# 
# names(datExpr)[moduleColors=="green"]
# 
# 
# #  Code chunk 9
# 
# # Create the starting data frame
# geneInfo0 = data.frame(gene = names(datExpr),
#                        #geneSymbol = annot$gene_symbol[probes2annot],
#                        #LocusLinkID = annot$LocusLinkID[probes2annot],
#                        moduleColor = moduleColors,
#                        geneTraitSignificance,
#                        GSPvalue)
# # Order modules by their significance for weight
# modOrder = order(-abs(cor(MEs, yan, use = "p")));
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
#   oldNames = names(geneInfo0)
#   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
#                          MMPvalue[, modOrder[mod]]);
#   names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.yan));
# geneInfo = geneInfo0[geneOrder, ]
# 
# 
# #  Code chunk 10
# 
# write.csv(geneInfo, file = "sandbox/wgcna/geneInfo-day1.csv")
