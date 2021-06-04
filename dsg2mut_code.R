#DSG2Mut Experiment
#Collaboration with Les Tung Lab
#Primary code author: Suraj Kannan
#Primary project lead: Robert Hawthorne
#May 31, 2021
#Document: Code + Figures 1.0

library(ggplot2)
library(pheatmap)
library(singleCellNet)
library(Matrix)
library(DESeq2)
library(pathview)

#Set the working directory
setwd("~/Documents/Research/RobHawthorne/")

#Load and clean data
data = read.table("~/Documents/Research/RobHawthorne/counts_table.txt", as.is = TRUE, row.names = 1, header = TRUE)
data = data[, 6:15]

#The featureCounts mapping was done to ensembl, but I generally prefer to work in gene ids (I know - scandalous!). There are many ways to do this conversion (for example, biomaRt); however, I have previously written a couple of my own custom scripts to do this. To get this part running, you'll need to grab the clean_nodatasets workspace - you can get the instructions to get this from my github, e.g. https://github.com/skannan4/cm-entropy-score.
load("~/Documents/Research/Reference/Cleanup/clean_nodatasets_060720.RData")
clean_data = rename_genes(data, species = "human")
colnames(clean_data) = substr(colnames(clean_data), 7, 10)

#Get the phenotype table
pheno = read.table("~/Documents/Research/RobHawthorne/phenotype_table.txt", as.is = TRUE, row.names = 1, header = TRUE)
pheno$depth = colSums(clean_data)
pheno$genes = colSums(clean_data > 0)
pheno$mito = mito(clean_data, species = "human")
clean_data_int = apply(clean_data, 2, as.integer) #DESeq likes to work with integers. The decimals in the counts table, by the way, come from the multimapping settings 
rownames(clean_data_int) = rownames(clean_data)

#Set up the DESeq object and run
dds = DESeqDataSetFromMatrix(countData = clean_data_int, colData = pheno, design = ~Nickname)
dds = DESeq(dds)

wt_vs_mutant = as.data.frame(results(dds, contrast=c("Nickname","wildtype","DSG2_mutant")))
wt_vs_mutant = wt_vs_mutant[!is.na(wt_vs_mutant$padj), ]
wt_vs_mutant = wt_vs_mutant[order(wt_vs_mutant$padj), ]
wt_vs_mutant_genes = rownames(wt_vs_mutant)[wt_vs_mutant$padj < 0.05]

#Some basic exploratory plots
rlog = rlogTransformation(dds)
plotPCA(rlog, intgroup = c("Nickname"))
pca = prcomp(t(assay(rlog)))
pcaData = as.data.frame(pca$x)
pcaData$Nickname = pheno$Nickname
pheatmap(assay(rlog)[wt_vs_mutant_genes, ], show_rownames = FALSE, scale = "row", annotation = pheno[, c(4, 6)])

#Making pathview figures - example here is for Nfkb pathway, but any pathway of interest can be used
pv.out <- pathview(gene.data = wt_vs_mutant[wt_vs_mutant$padj < 0.05, 2, drop = FALSE], pathway.id = "04064", species = "hsa", out.suffix = "nfkb_wt_mutant", gene.idtype = "symbol", low = list(gene = "red", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "green", cpd = "yellow"))

###Making a pretty heatmap
breaks_heatmap = unique(c(seq(-3, 0, 0.1), 0, seq(0, 3, 0.1)))
pheno$Condition = c(rep("Control", 5), rep("DSG2mut", 5))
temp_cols = c("#3C4A80", "#DC3B24")
names(temp_cols) = c("Control", "DSG2mut")
colors_heatmap = c(colorRampPalette(c("blue", "white"))(30), "white", colorRampPalette(c("white", "red"))(30))
goi = read.table("goi.txt", as.is = TRUE)$V1
rlog = rlog[rowSds(assay(rlog)) != 0, ]
pheatmap(assay(rlog)[goi[goi %in% wt_vs_mutant_genes], pheno$Nickname %in% c("wildtype", "DSG2_mutant")], show_rownames = TRUE, show_colnames = FALSE, scale = "row", annotation = pheno[pheno$Nickname %in% c("wildtype", "DSG2_mutant"), c(5), drop = FALSE], color = colors_heatmap, breaks = breaks_heatmap, cellheight = 10, cellwidth = 10, cluster_cols = FALSE, annotation_colors = list(Condition = temp_cols))