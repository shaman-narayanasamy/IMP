#!/bin/R

###################################################################################################
## Load required packages
###################################################################################################

require(genomeIntervals)

require(checkpoint)
checkpoint('2015-04-27', scanForPackages=FALSE)

require(ggplot2)
require(gtools)
require(data.table)
require(reshape)
require(grid)
require(grDevices)
require(stringr)
require(xtable)
require(beanplot)
require(psych)

###################################################################################################
## Read in arguments
###################################################################################################
# User should input the three different files
# Uncomment when ready

print("START: Reading arguments")

args		    <- commandArgs(trailingOnly = TRUE)
out_dir		    <- args[1] # Output directory
MG.read.count_file  <- args[2]
MT.read.count_file  <- args[3]
MG.map.summary_file <- args[4]
MT.map.summary_file <- args[5]
MG.cov_file	    <- args[6]
MT.cov_file	    <- args[7]
MG.depth_file	    <- args[8]
MT.depth_file	    <- args[9]
MG.var_file	    <- args[10]
MT.var_file	    <- args[11]
GC.dat_file	    <- args[12]
coords_file	    <- args[13]
annot_file	    <- args[14]
nucmer_file	    <- args[15]

print("DONE: Reading arguments")

###################################################################################################
## Initialize functions for various calculations and normalizations
###################################################################################################
print("START: Reading functions")
source("src/IMP_plot_functions.R")
print("DONE: Reading functions")

###################################################################################################
## Read in the necessary input files
###################################################################################################
print("START: Reading data")

## Read counts MG
print("Read in MG read count file")
MG.read.count <- read.table(MG.read.count_file)
MG.read.count <- as.data.frame(cbind(as.character(MG.read.count$V2)[-nrow(MG.read.count)], MG.read.count$V1[-nrow(MG.read.count)]))
colnames(MG.read.count) <- c("file", "count")
MG.read.count$file <- as.character(MG.read.count$file)
# Divide by four for total sequences because original file is a line count
MG.read.count$count <- as.numeric(as.character(MG.read.count$count))/4

## Read counts MT
print("Read in MT read count file")
MT.read.count <- read.table(MT.read.count_file)
MT.read.count <- as.data.frame(cbind(as.character(MT.read.count$V2)[-nrow(MT.read.count)], MT.read.count$V1[-nrow(MT.read.count)]))
colnames(MT.read.count) <- c("file", "count")
MT.read.count$file <- as.character(MT.read.count$file)
# Divide by four for total sequences because original file is a line count
MT.read.count$count <- as.numeric(as.character(MT.read.count$count))/4

## Read MG mapping summary file (samtools flagstat)
print("Read in MG mapping summary file")
MG.map.summary <- read.table(MG.map.summary_file)
MG.map.summary <- cbind(
		     c("Total",
		       "Duplicates",
		       "Mapped",
		       "Paired",
		       "Read1",
		       "Read2",
		       "Properly paired",
		       "With itself and mate",
		       "Singletons",
		       "Mate mapped to different contig",
		       "Mate mapped to different contig mapQ>=5"
		       ),
		     MG.map.summary
		     )
colnames(MG.map.summary) <- c("Mapping", "Reads")

## Read MT mapping summary file (samtools flagstat)
print("Read in MT mapping summary file")
MT.map.summary <- read.table(MT.map.summary_file)
MT.map.summary <- cbind(
		     c("Total",
		       "Duplicates",
		       "Mapped",
		       "Paired",
		       "Read1",
		       "Read2",
		       "Properly paired",
		       "With itself and mate",
		       "Singletons",
		       "Mate mapped to different contig",
		       "Mate mapped to different contig mapQ>=5"
		       ),
		     MT.map.summary
		     )
colnames(MT.map.summary) <- c("Mapping", "Reads")

## MG coverage (samtools idxstat output)
print("Read in MG coverage file")
MG.cov <- read.table(MG.cov_file, colClasses=c("factor", rep("numeric", 6 )))[,c(1,4:7)]
colnames(MG.cov) <- c("contig", "MG_reads", "MG_bases_covered", "length", "MG_cov")

## MT coverage (samtools idxstat output)
print("Read in MT coverage file")
MT.cov <- read.table(MT.cov_file, colClasses=c("factor", rep("numeric", 6 )))[,c(1,4:7)]
colnames(MT.cov) <- c("contig", "MT_reads", "MT_bases_covered", "length", "MT_cov")

## MG depth (bedtools coverageBed output)
print("Read in MG depth file")
MG.depth <- read.table(MG.depth_file, colClasses=c("factor", "numeric"))
colnames(MG.depth) <- c("contig", "MG_depth")

## MT depth (bedtools coverageBed output)
print("Read in MT depth file")
MT.depth <- read.table(MT.depth_file, colClasses=c("factor", "numeric"))
colnames(MT.depth) <- c("contig", "MT_depth")

## MG variant calls
print("Read in MG variant calls")
MG.var <- as.data.frame(table(read.table(MG.var_file)[,1]))
colnames(MG.var) <- c("contig", "MG_var")

## MT variant calls
print("Read in MT variant calls")
MT.var <- as.data.frame(table(read.table(MT.var_file)[,1]))
colnames(MT.var) <- c("contig", "MT_var")

## GC content file
print("Read in GC percentage file")
GC.dat <- read.delim(GC.dat_file)
colnames(GC.dat) <- c("contig", "GC")

## gff annotation file
print("Read in gff3 annotation file")
annot <- readGff3(annot_file)

# extract only interesting information
print("Processing gff3 annotation file")
annot.1 <- as.data.frame(cbind(as.character(annot@annotation$seq_name), annot@.Data,
      str_split_fixed(str_split_fixed(getGffAttribute(annot, "inference"), ",", 2)[,2], ":", 3)[,2]))
colnames(annot.1) <- c("contig", "start", "end", "database")
annot.1$database[annot.1$database==""]=NA

annot.1[,2:3] <- sapply(annot.1[2:3], as.character)
annot.1[,2:3] <- sapply(annot.1[2:3], as.numeric)

# calculate no. of gene length
print("Calculating no. of genes per contig")
annot.2 <- as.data.frame(table(annot.1$contig))
colnames(annot.2)[1:ncol(annot.2)] <- c("contig","no_of_genes")

# calculate gene lengths
print("Calculating length of genes from gff3 file")
annot.1 <- as.data.frame(cbind(annot.1, annot.1$end - annot.1$start + 1))
colnames(annot.1)[ncol(annot.1)] <- "gene_length"

# aggregate table and calculate total gene lengths within contig
print("Calculating coding density of contigs")
annot.2 <- cbind(annot.2, aggregate(gene_length~contig, data=annot.1, FUN=sum)[,2])
colnames(annot.2)[ncol(annot.2)] <- "total_gene_length"

# create annotation table
print("Creating annotation table")
annot.3 <- as.data.frame.matrix(table(annot.1[,c(1,4)]))[,-1]
annot.3 <- cbind(rownames(annot.3), annot.3, rowSums(annot.3[,c(2:ncol(annot.3))]))
rownames(annot.3) <- NULL
colnames(annot.3)[c(1, ncol(annot.3))] <- c("contig", "all_annotations")

annot_final <- merge(annot.2, annot.3, by="contig")

# vizbin points, the contig names are usually not included in this file,
# concatenate it before reading in
print("Reading in vizbin coordinates")
coords <- read.table(coords_file, colClasses=c("factor", "numeric", "numeric"),
		     sep="\t", col.names=c("contig", "x", "y"))

print("Reading in nucmer results")
nucmer_res <- read.table(nucmer_file, header=F)
colnames(nucmer_res) <- c("ref_start", "ref_end", "query_start", "query_end", "ref_align_len",
			  "query_align_len", "identity", "ref_id", "contig")

print("DONE: Reading data")
###################################################################################################
## Merge the data sets without the vizbin coordinates
###################################################################################################
print("START: Merging data")

all.dat <- merge(MG.cov, MT.cov, by=c("contig", "length"), all=T)
all.dat <- merge(all.dat, GC.dat, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MG.depth, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MT.depth, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MG.var, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MT.var, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, annot_final, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, nucmer_res, by=c("contig"), all=T, incomparables=NA)

print("DONE: Merging data")
###################################################################################################
## Perform calculations and append it to the full table
###################################################################################################
# Calculate and merge data
print("START: Performing calculations")

# Compute gene density
print("Computing gene density")
all.dat$gene_dens <- gene_density(all.dat$no_of_genes, all.dat$length)

# Compute coding density
print("Computing coding density")
all.dat$coding_dens <- coding_density(all.dat$total_gene_length, all.dat$length)

# Compute RPKM values
print("Computing MG RPKM values")
all.dat$MG_rpkm <- contig_rpkm(all.dat$MG_reads, all.dat$length)
print("Computing MT RPKM values")
all.dat$MT_rpkm <- contig_rpkm(all.dat$MT_reads, all.dat$length)

# Compute variant densities
print("Computing MG variant densities")
all.dat$MG_var_dens <-  var_density(all.dat$MG_var, all.dat$length, all.dat$MG_reads)
print("Computing MT variant densities")
all.dat$MT_var_dens <-  var_density(all.dat$MT_var, all.dat$length, all.dat$MT_reads)

# Compute ratio values
print("Computing coverage ratio")
all.dat$cov_ratio <-  all.dat$MT_cov / all.dat$MG_cov
print("Computing depth ratio")
all.dat$depth_ratio <-  all.dat$MT_depth / all.dat$MG_depth
print("Computing RPKM ratio")
all.dat$rpkm_ratio <- all.dat$MT_rpkm / all.dat$MG_rpkm
print("Computing variant ratio")
all.dat$var_ratio <- all.dat$MT_var_dens / all.dat$MG_var_dens

# Compute log transformation of ratio values
print("Log transforming some of the data")
all.dat$log_cov_ratio <- log10(all.dat$cov_ratio)
all.dat$log_depth_ratio <- log10(all.dat$depth_ratio)
all.dat$log_rpkm_ratio <- log10(all.dat$rpkm_ratio)
all.dat$log_var_ratio <- log10(all.dat$var_ratio)

# Compute query coverage against reference genomes
print("Computing query coverage of contigs against reference genomes")
all.dat$query_cov <- all.dat$query_align_len / all.dat$length * 100
all.dat$query_cov[na.omit(all.dat$query_cov > 100)] = 100

print("DONE: Performing calculations")
###################################################################################################
## Organize filtering statistics and create table
###################################################################################################
print("START: Organizing preprocessing stats")

## Organize MG filtering procedures
print("Printing metagenomic filtering statistics")
# Obtain file names without path
MG.read.count$file <- unlist(lapply(MG.read.count$file, get_file_name))

# Obtain string length
MG.read.count <- cbind(MG.read.count,
		       unlist(lapply(MG.read.count$file, nchar)))
colnames(MG.read.count)[ncol(MG.read.count)] <- "string_len"

# Order the file names according to the string length
MG.read.count <- MG.read.count[order(MG.read.count$string_len),]

# Create additional column for fastq file types (paired or singletons)
MG.read.count <- cbind(MG.read.count, rep(NA, nrow(MG.read.count)))
colnames(MG.read.count)[ncol(MG.read.count)] <- "type"
MG.read.count$type[grep("R[12]", MG.read.count$file)] <- "paired-end"
MG.read.count$type[grep("SE", MG.read.count$file)] <- "singletons"

# Create new columns with filtering type (within file names)
MG.read.count <- cbind(MG.read.count,
		      unlist(lapply(MG.read.count$file, filtering)))
colnames(MG.read.count)[ncol(MG.read.count)] <- "filtering"
MG.read.count$filtering <- as.character(MG.read.count$filtering)

# Rename filtering steps
MG.read.count$filtering[ # Unfiltered/raw
			grep("R[12]", MG.read.count$filtering)] <-
			   "unfiltered (raw)"
MG.read.count$filtering <- # Remove _filt extension in file name
   gsub("_filt", "", MG.read.count$filtering, ignore.case=TRUE)
MG.read.count$filtering <- # Remove _filt extension in file name
   gsub("trimmed", "adapter trimmed and quality filtered", MG.read.count$filtering, ignore.case=TRUE)
MG.read.count$filtering <- # Convert uniq to unique
   gsub("uniq", "collapsed unique", MG.read.count$filtering, ignore.case=TRUE)
MG.read.count$filtering <- # Convert hg to "human sequences"
    gsub("hg\\d+", "human sequences filtered", MG.read.count$filtering,
	ignore.case=TRUE)

# Summarize table for final report
MG.read.count.final <- unique(MG.read.count[,c(5,4,2)])
MG.read.count.final$count <- as.integer(MG.read.count.final$count)

# Write out table in html and tab separated files
sink(name_plot("MG.read_stats.html"))
print(xtable(MG.read.count.final, html.table.attributes=""), type = "html")
sink()

write.table(MG.read.count.final, name_plot("MG.read_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

####################################################################
## Organize MT filtering procedures
####################################################################

print("Printing metatranscriptomic filtering statistics")
# Obtain file names without path
MT.read.count$file <- unlist(lapply(MT.read.count$file, get_file_name))

# Obtain string length
MT.read.count <- cbind(MT.read.count,
		       unlist(lapply(MT.read.count$file, nchar)))
colnames(MT.read.count)[ncol(MT.read.count)] <- "string_len"

# Order the file names according to the string length
MT.read.count <- MT.read.count[order(MT.read.count$string_len),]

# Create additional column for fastq file types (paired or singletons)
MT.read.count <- cbind(MT.read.count, rep(NA, nrow(MT.read.count)))
colnames(MT.read.count)[ncol(MT.read.count)] <- "type"
MT.read.count$type[grep("R[12]", MT.read.count$file)] <- "paired-end"
MT.read.count$type[grep("SE", MT.read.count$file)] <- "singletons"

# Create new columns with filtering type (within file names)
MT.read.count <- cbind(MT.read.count,
		      unlist(lapply(MT.read.count$file, filtering)))
colnames(MT.read.count)[ncol(MT.read.count)] <- "filtering"
MT.read.count$filtering <- as.character(MT.read.count$filtering)

# Rename filtering steps
MT.read.count$filtering[ # Unfiltered/raw
			grep("R[12]", MT.read.count$filtering)] <-
			   "unfiltered (raw)"
MT.read.count$filtering <- # Remove _filt extension in file name
   gsub("_filt", "", MT.read.count$filtering, ignore.case=TRUE)
MT.read.count$filtering <- # Remove _filt extension in file name
   gsub("trimmed", "adapter trimmed and quality filtered", MT.read.count$filtering, ignore.case=TRUE)
MT.read.count$filtering <- # Convert uniq to unique
   gsub("uniq", "collapsed unique", MT.read.count$filtering, ignore.case=TRUE)
MT.read.count$filtering <- # Convert rna to "rrna sequences"
    gsub("rna", "rrna sequences filtered", MT.read.count$filtering,
	ignore.case=TRUE)
MT.read.count$filtering <- # Convert hg to "human sequences"
    gsub("hg\\d+", "human sequences filtered", MT.read.count$filtering,
	ignore.case=TRUE)

# Summarize table for final report
MT.read.count.final <- unique(MT.read.count[,c(5,4,2)])
MT.read.count.final$count <- as.integer(MT.read.count.final$count)

# Write out table in html and tab separated files
sink(name_plot("MT.read_stats.html"))
print(xtable(MT.read.count.final, html.table.attributes=""), type = "html")
sink()

write.table(MT.read.count.final, name_plot("MT.read_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

print("DONE: Organizing preprocessing stats")

###################################################################################################
## Incorporate coordinates from vizbin
###################################################################################################
print("START: Incorporating VizBin data")

vb_dat <- merge(all.dat, coords, by=c("contig"), all=F, incomparables=NA)
vb_dat <- vb_dat[!is.na(vb_dat$x),]

print("Processing VizBin data")
# Handle missing and infinite values
vb_dat[is.na(vb_dat)] <- 0

# Remove outliers and infinite values
print("Handling outliers")

vb_dat$MG_var_dens <- outliers(vb_dat$MG_var_dens,2)
vb_dat$MT_var_dens <- outliers(vb_dat$MT_var_dens,2)
vb_dat$MG_depth <- outliers(vb_dat$MG_depth,2)
vb_dat$MT_depth <- outliers(vb_dat$MT_depth,2)
vb_dat$MG_rpkm <- outliers(vb_dat$MG_rpkm,2)
vb_dat$MT_rpkm <- outliers(vb_dat$MT_rpkm,2)
vb_dat$cov_ratio <- outliers(vb_dat$cov_ratio,2)
vb_dat$depth_ratio <- outliers(vb_dat$depth_ratio,2)
vb_dat$rpkm_ratio <- outliers(vb_dat$rpkm_ratio,2)
vb_dat$var_ratio <- outliers(vb_dat$var_ratio,2)
vb_dat$log_cov_ratio <- outliers(vb_dat$log_cov_ratio,2)
vb_dat$log_depth_ratio <- outliers(vb_dat$log_depth_ratio,2)
vb_dat$log_rpkm_ratio <- outliers(vb_dat$log_rpkm_ratio,2)
vb_dat$log_var_ratio <- outliers(vb_dat$log_var_ratio,2)

print("DONE: Incorporating VizBin data")
####################################################################
## ASSEMBLY STATISTICS AND VISUALIZATIONS
####################################################################
print("START: Visualizing")
## Output filtering statistics table (text and tsv)

## Output mapping stats table
print("Print metagenomic mapping statistics table")
sink(name_plot("MG_mapping_stats.html"))
print(xtable(MG.map.summary, html.table.attributes=""), type="html")
sink()

write.table(MG.map.summary, name_plot("MG_mapping_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

## Output mapping stats table
print("Print metatranscriptomics mapping statistics table")
sink(name_plot("MT_mapping_stats.html"))
print(xtable(MT.map.summary, html.table.attributes=""), type="html")
sink()

write.table(MT.map.summary, name_plot("MT_mapping_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

## Plot standard vizbin scatter plot (non included)
print("Generating standard vizbin plot")
png(name_plot("IMP-vizbin_standard.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point() +
ggtitle("Standard VizBin") +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length info
print("Generating vizbin plot with contig length information")
png(name_plot("IMP-vizbin_length.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length))) +
guides(size=guide_legend(title=log10len)) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and GC content
print("Generating vizbin plot for GC content")
png(name_plot("IMP-vizbin_length_GC.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length), colour=GC, order=GC),alpha=0.5) +
scale_colour_gradientn(colours=rainbow(7), guide="colourbar", guide_legend(title="% G+C")) +
guides(size=guide_legend(title=log10len)) +
theme_nothing()
dev.off()

## Vizbin plot for quast results
print("Generating vizbin plot for quast results")
nref <- length(unique(vb_dat$ref_id))
png(name_plot("IMP-vizbin_quast.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=query_cov,
	       colour=ref_id,
	       order=length,
	       alpha=identity)) +
		   guides(size=guide_legend(title="Query coverage"),
		   colour=guide_legend(title="Reference"),
		   alpha=guide_legend(title="% identity")
		   )+
theme_gray()
dev.off()

####################################################################
## MAPPING STATISTICS AND VISUALIZATIONS
####################################################################

## Plot mappable reads density
print("Generating mapped reads plot")
var1 <-log10(c(all.dat$MG_reads,all.dat$MT_reads))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MT",nrow(all.dat)))
mapped_reads<-data.frame(var1,var2)

png(name_plot("IMP-reads_density.png"), width=700, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= mapped_reads,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Mappable reads", ylab=expression(log[10]*~"count"))
legend("bottomleft", fill =c("blue", "red"), legend = c("MG", "MT"))
dev.off()

## Plot rpkm density
print("Generating rpkm plot")
var1 <-log10(c(all.dat$MG_rpkm,all.dat$MT_rpkm))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MT",nrow(all.dat)))
rpkm<-data.frame(var1,var2)

png(name_plot("IMP-rpkm_density.png"), width=700, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= rpkm,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="RPKM", ylab=expression(log[10]*~"RPKM"))
legend("bottomleft", fill =c("blue", "red"), legend = c("MG", "MT"))
dev.off()

## Plot coverage density
print("Generating coverage plot")
var1 <-c(all.dat$MG_cov,all.dat$MT_cov)
var2 <- c(rep("MG",nrow(all.dat)),rep("MT",nrow(all.dat)))
coverage<-data.frame(var1,var2)

png(name_plot("IMP-coverage_density.png"), width=700, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= coverage,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Coverage", ylab="fraction")
legend("bottomleft", fill =c("blue", "red"), legend = c("MG", "MT"))
dev.off()

## Plot depth density
print("Generating depth density plot")
var1 <-log10(c(all.dat$MG_depth,all.dat$MT_depth))
var2 <- c(rep("MG",nrow(all.dat)),rep("MT",nrow(all.dat)))
depth<-data.frame(var1,var2)

png(name_plot("IMP-depth_density.png"), width=700, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= depth,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Depth", ylab=expression(log[10]*~"avg. depth"))
legend("bottomleft", fill =c("blue", "red"), legend = c("MG", "MT"))
dev.off()

## Plot vizbin scatter plot with length and MG coverage info
print("Generating vizbin plot for metagenomic coverage")
png(name_plot("IMP-vizbin_length_MGcov.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MG_cov, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "coverage"))))
      ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MT coverage info
print("Generating vizbin plot for metatranscriptomic coverage")
png(name_plot("IMP-vizbin_length_MTcov.png"),width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="red", aes(alpha=MT_cov, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metatranscriptomic", "coverage"))))
      ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MG depth info
print("Generating vizbin plot for metagenomic depth")
png(name_plot("IMP-vizbin_length_MGdepth.png"),width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MG_depth, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "depth"))))
      ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MT depth info
print("Generating vizbin plot for metatranscriptomic depth")
png(name_plot("IMP-vizbin_length_MTdepth.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="red", aes(alpha=MT_depth, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metatranscriptomic", "depth"))))
      ) +
theme_nothing()
dev.off()

####################################################################
## VARIANT STATISTICS AND VISUALIZATIONS
####################################################################
## Plot variant count
var1 <-log10(c(all.dat$MG_var,all.dat$MT_var))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MT",nrow(all.dat)))
variant_count <- data.frame(var1,var2)

print("Generating variant count plots")
png(name_plot("IMP-var_count.png") ,width=700, height=700)

par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= variant_count,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
main="Variant count", ylab=expression(log[10]*~count), bw="nrd")
legend("bottomleft", fill =c("blue", "red"), legend = c("MG", "MT"))
dev.off()

## Plot variant density
var1 <-c(all.dat$MG_var_dens,all.dat$MT_var_dens)
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MT",nrow(all.dat)))
variant_density<-data.frame(var1,var2)

print("Generating variant density plots")
png(name_plot("IMP-var_density.png") ,width=700, height=700)

par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= variant_density,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
main="Variant density", ylab="count / RPKM", bw="nrd")
legend("bottomleft", fill =c("blue", "red"), legend = c("MG", "MT"))
dev.off()

## Plot vizbin scatter plot with length and MG variant info
# Create label
MG_var_label <- expression(bold(frac(variants[MG]/kb, MG[rpkm])))

print("Generating vizbin plot for metagenomic variant density")
png(name_plot("IMP-vizbin_length_MGvardens.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(colour=MG_var_dens, size=log10(length), order=MG_var_dens), alpha=0.75) +
scale_colour_gradient(high="black", low="cyan") +
       guides(size=guide_legend(title=log10len),
	      colour=guide_colourbar(title=MG_var_label)
       ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MT variant info
# Create label
MT_var_label <- expression(bold(frac(variants[MT]/kb, MT[rpkm])))

print("Generating vizbin plot for metatranscriptomic variant density")
png(name_plot("IMP-vizbin_length_MTvardens.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(colour=MT_var_dens, size=log10(length), order=MT_var_dens), alpha=0.75) +
scale_colour_gradient(high="black", low="magenta") +
       guides(size=guide_legend(title=log10len),
	      colour=guide_colourbar(title=MT_var_label)
       ) +
theme_nothing()
dev.off()

####################################################################
## ANNOTATION STATISTICS AND VISUALIZATIONS
####################################################################
## Vizbin plot with length and raw number of genes

print("Generating vizbin plot for number of genes")
png(name_plot("IMP-vizbin_length_geneCount.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(colour=no_of_genes, size=log10(length), order=no_of_genes), alpha=0.75) +
scale_colour_gradientn(colours=topo.colors(max(vb_dat$no_of_genes)),
		      guide="colourbar",
		      guide_legend(title="Gene count")
		      ) +
       guides(size=guide_legend(title=log10len)
       ) +
theme_nothing()
dev.off()

## Vizbin plot with length and gene density
print("Generating vizbin plot for gene density")
png(name_plot("IMP-vizbin_length_geneDensity.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(colour=gene_dens, size=log10(length), order=gene_dens), alpha=0.75) +
scale_colour_gradientn(colours=topo.colors(500),
		      guide="colourbar",
		      guide_legend(title="Gene density\n (genes/1kb)")
		      ) +
       guides(size=guide_legend(title=log10len)
       ) +
theme_nothing()
dev.off()

####################################################################
## RATIO STATISTICS AND VISUALIZATION
####################################################################
print("Generate statistics and visualizations")

## Vizbin plot for coverage ratio
print("Generating vizbin plot for coverage ratios")
png(name_plot("IMP-vizbin_length_covRatio.png"), width=700, height=700)
covRatio_label <- expression(bold(paste(log[10], frac(cov[MT], cov[MG]))))
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length),
	       colour=log_cov_ratio,
	       order=log_cov_ratio), alpha=0.5) +
scale_colour_gradientn(colours=rev(heat.colors(100))) +
		       guides(size=guide_legend(title=log10len),
		       colour=guide_colourbar(title=covRatio_label)
		       )+
theme_black()
dev.off()

## Vizbin plot for depth ratio
print("Generating vizbin plot for depth ratios")
png(name_plot("IMP-vizbin_length_depthRatio.png"), width=700, height=700)
depthRatio_label <- expression(bold(paste(log[10], frac(depth[MT], depth[MG]))))
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length),
	       colour=log_depth_ratio,
	       order=log_depth_ratio), alpha=0.5) +
scale_colour_gradientn(colours=rev(heat.colors(1000))) +
		   guides(size=guide_legend(title=log10len),
		   colour=guide_colourbar(title=depthRatio_label)
		   )+
theme_black()
dev.off()

## Vizbin plot for rpkm ratio
print("Generating vizbin plot for rpkm ratios")
png(name_plot("IMP-vizbin_length_rpkmRatio.png"), width=700, height=700)
rpkmRatio_label <- expression(bold(paste(log[10], frac(rpkm[MT], rpkm[MG]))))
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length),
	       colour=log_rpkm_ratio,
	       order=log_rpkm_ratio),
	   alpha=0.5) +
scale_colour_gradientn(colours=rev(heat.colors(1000))) +
		   guides(size=guide_legend(title=log10len),
		   colour=guide_colourbar(title=rpkmRatio_label)
		   )+
theme_black()
dev.off()

print("DONE: Visualizing")

####################################################################
## Save the R workspace
####################################################################

print("START: Saving R image: MGMT_results.Rdat")
save.image(name_plot("MGMT_results.Rdat"))
print("DONE: Saving R image: MGMT_results.Rdat")

