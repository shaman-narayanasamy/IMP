#!/bin/R

### Arguments
args <- commandArgs(TRUE)
WKSPC <- args[1] # Path to workspace
BED <- args[2] # Path to links
PID <- args[3] # Essential gene hits
SAMPLE <- args[4] # Sample number

# Load workspace
load(WKSPC)

# List essential genes
markerGenes <- c("TIGR01011", "TIGR01049", 
		 "TIGR01169", "TIGR00487", 
		 "TIGR01044", "TIGR00959", 
		 "PF00573.18", "TIGR01171", 
		 "PF00380.15", "PF00297.18", 
		 "TIGR00472", "TIGR01067", 
		 "TIGR01021", "TIGR01050",
		 "TIGR01029", "TIGR01164",
		 "PF00416.18", "TIGR00468",
		 "TIGR01071", "PF00276.16",
		 "PF00347.19", "TIGR01632",
		 "PF00281.15", "TIGR00981",
		 "TIGR00012", "TIGR01009",
		 "PF00411.15", "PF00466.16",
		 "PF00410.15", "TIGR00060",
		 "TIGR00952", "PF00366.16",
		 "TIGR01066", "TIGR01079",
		 "TIGR02013")

# Read essential gene links
essLinks <- read.table(BED, header=F)[,c(1,4)]
colnames(essLinks) <- c("contig", "gene_ID")

# Read essential gene annotation
essAnnot <- read.table(PID, header=F)[,c(1,3,4)]
colnames(essAnnot) <- c("gene_ID", "gene_function", "gene_name")

## Merge essential gene results with contigs 
essAll <- merge(essAnnot, essLinks, by="gene_ID", all=F) 
 
## Merge essential gene results with clusters 
clusterRes <- merge(clusterRes, essAll, by="contig", all=F)

clusterRes <- cbind(sample=rep(SAMPLE, nrow(clusterRes)), clusterRes)

## Rearrange the table in following order of columns: gene ID, sample ID, cluster ID, gene name.
clusterRes <- clusterRes[,c("gene_ID", "sample", "cluster", "gene_name")]

## Retain only marker genes
clusterRes <- clusterRes[clusterRes$gene_name%in%markerGenes,]

## Remove N and B clusters
clusterRes <- clusterRes[-c(grep("^[NB]", clusterRes$cluster, perl = TRUE)),]
clusterRes <- unique(clusterRes)

## Write out table
write.table(clusterRes, "Binning/essMarkerGenes/markersAll.tsv", quote=F, row.names=F, sep="\t", col.names=F)

## Write out separate tables for the cluster of genes
invisible(
lapply(markerGenes, 
       FUN=function(x){write.table(subset(clusterRes, gene_name==x, drop=FALSE),
                                   paste("Binning/essMarkerGenes/marker-", x, ".tsv", sep=""),
                                   quote=F, row.names=F, sep="\t", col.names=F 
                                   )}
       )
)
