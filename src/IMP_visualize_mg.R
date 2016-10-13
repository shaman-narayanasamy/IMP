#!/bin/R

###################################################################################################
## Load required packages
###################################################################################################

print("Loading required R libraries")
require(genomeIntervals)

require(checkpoint)
checkpoint('2016-06-20', scanForPackages=FALSE, checkpointLocation="~/lib", project="~/lib")

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
MG.map.summary_file <- args[3]
MG.cov_file	    <- args[4]
MG.depth_file	    <- args[5]
MG.var_file	    <- args[6]
GC.dat_file	    <- args[7]
annot_file	    <- args[8]
coords_file	    <- args[9]
nucmer_file	    <- args[10]
function_script	    <- args[11]
print("DONE: Reading arguments")

###################################################################################################
## Initialize functions for various calculations and normalizations
###################################################################################################

print("START: Reading functions")
source(function_script)
print("DONE: Reading functions")

###################################################################################################
## Read in the necessary input files
###################################################################################################

print("Reading input files...")

## Read counts MG
print("Read in MG read count file")
MG.read.count <- read.table(MG.read.count_file)
MG.read.count <- as.data.frame(cbind(as.character(MG.read.count$V2)[-nrow(MG.read.count)], MG.read.count$V1[-nrow(MG.read.count)]))
colnames(MG.read.count) <- c("file", "count")
MG.read.count$file <- as.character(MG.read.count$file)
# Divide by four for total sequences because original file is a line count
MG.read.count$count <- as.numeric(as.character(MG.read.count$count))/4

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

## MG coverage (samtools idxstat output)
print("Read in MG coverage file")
MG.cov <- read.table(MG.cov_file, colClasses=c("factor", rep("numeric", 6 )))[,c(1,4:7)]
colnames(MG.cov) <- c("contig", "MG_reads", "MG_bases_covered", "length", "MG_cov")

## MG depth (bedtools coverageBed output)
print("Read in MG depth file")
MG.depth <- read.table(MG.depth_file, colClasses=c("factor", "numeric"))
colnames(MG.depth) <- c("contig", "MG_depth")

## MG variant calls
print("Read in MG variant calls")
MG.var <- as.data.frame(table(read.table(MG.var_file)[,1]))
colnames(MG.var) <- c("contig", "MG_var")

## GC content file
print("Read in GC percentage file")
GC.dat <- read.delim(GC.dat_file)
colnames(GC.dat) <- c("contig", "GC")

## gff annotation file
print("Read in gff3 annotation file")
annot <- readGff3(annot_file, isRightOpen=T)

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
# Create temporary table
total_gene_length <- aggregate(gene_length~contig, data=annot.1, FUN=sum)
colnames(total_gene_length)[ncol(total_gene_length)] <- "total_gene_length"
annot.2 <- merge(annot.2, 
		 total_gene_length,
		 by="contig", 
		 all.y=FALSE)

# remove temporary table
rm(total_gene_length)

# create annotation table
print("Creating annotation table")
annot.3 <- as.data.frame.matrix(table(annot.1[,c(1,4)]))[,-1]
annot.3 <- cbind(rownames(annot.3), annot.3, rowSums(annot.3[,c(1:ncol(annot.3))]))
rownames(annot.3) <- NULL
colnames(annot.3)[c(1, ncol(annot.3))] <- c("contig", "all_annotations")


if (is.null(ncol(annot.3)))
{
    print("Skipping annotation count")
    annot.4 <- annot.2
}else{
    annot.3 <- cbind(annot.3, rowSums(annot.3[,c(3:ncol(annot.3))]))
    colnames(annot.3)[c(1, ncol(annot.3))] <- c("contig", "all_annotations")
    annot.4 <- merge(annot.2, annot.3, by="contig")
}


# vizbin points, the contig names are usually not included in this file,
# concatenate it before reading in
print("Reading in vizbin coordinates")
coords <- read.table(coords_file, colClasses=c("factor", "numeric", "numeric"),
		     sep="\t", col.names=c("contig", "x", "y"))

# Read in nucmer output from metaquast
print("Reading in nucmer results")
nucmer_try <- try(read.table(nucmer_file, header=F), silent=T)
if(inherits(nucmer_try, "try-error")){
  print("WARNING: Nucmer file empty. No taxanomy was assigned to contigs")
  nucmer_res <- read.table(text = "",
			   col.names = c("ref_start", "ref_end", "query_start", "query_end",
					 "ref_align_len", "query_align_len", "identity",
					 "ref_id", "contig")
	     )
}else{
  nucmer_res <- read.table(nucmer_file, header=F)
  colnames(nucmer_res) <- c("ref_start", "ref_end", "query_start", "query_end", "ref_align_len",
			  "query_align_len", "identity", "ref_id", "contig")
}

print("DONE: Reading data")
###################################################################################################
## Merge the data sets without the vizbin coordinates
###################################################################################################

print("START: Merging data")
all.dat <- merge(MG.cov, GC.dat, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MG.depth, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MG.var, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, annot.4, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, nucmer_res, by=c("contig"), all=T, incomparables=NA)
print("DONE: Merging data")

###################################################################################################
## Perform calculations and append it to the full table
###################################################################################################
# Calculate and merge data
print("Perform calculations")
save.image(name_plot("mg_results.Rdat"))
# Get new column numbers
newcols <- ncol(all.dat) + 1
all.dat <- cbind(all.dat,
      gene_density(all.dat$no_of_genes, all.dat$length),
      coding_density(all.dat$total_gene_length, all.dat$length),
      contig_rpkm(all.dat$MG_reads, all.dat$length),
      var_density(all.dat$MG_var, all.dat$length, all.dat$MG_reads)
)
colnames(all.dat)[newcols:ncol(all.dat)] <- c(
			      "gene_dens",
			      "coding_dens",
			      "MG_rpkm",
			      "MG_var_dens"
			      )
# Get new column numbers
newcols <- ncol(all.dat) + 1

###################################################################################################
## Organize filtering statistics and create table
###################################################################################################
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
sink(name_plot("mg.read_stats.html"))
print(xtable(MG.read.count.final, html.table.attributes=""), type = "html")
sink()

write.table(MG.read.count.final, name_plot("mg.read_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

###################################################################################################
## Calculate assembly statistics and create table
###################################################################################################
## Prepare table for assembly statistics
print("Printing assembly statistics")
assembly.stats <- cbind(c("All contigs", "Contigs >= 500", "Contigs >= 1000"),
		   rbind(
      get_stats(all.dat),
      get_stats(all.dat[all.dat$length>=500,]),
      get_stats(all.dat[all.dat$length>=1000,])
      ))

colnames(assembly.stats)<-c("Contig_set",
			"No_of_contigs",
			"N50",
			"Max_size",
			"Mean_size",
			"Median_size",
			"Total_length")

assembly.stats <- as.data.frame(assembly.stats)
assembly.stats[,2:6] <- sapply(assembly.stats[2:6], as.character)
assembly.stats[,2:6] <- sapply(assembly.stats[2:6], as.numeric)

###################################################################################################
## Incorporate coordinates from vizbin
###################################################################################################
print("Incorporating vizbin coordinates")
vb_dat <- merge(all.dat, coords, by=c("contig"), all=F, incomparables=NA)
vb_dat <- vb_dat[!is.na(vb_dat$x),]

#write.table(vb_dat, "final.contig.merged_min1000_info_raw.txt", sep="\t", quote=F, row.names=F)
# Handle missing and infinite values
vb_dat[is.na(vb_dat)] <- 0


print("Processing vizbin-based data")
# Remove outliers and infinite values
vb_dat$MG_var_dens <- outliers(vb_dat$MG_var_dens,2)
vb_dat$MG_depth <- outliers(vb_dat$MG_depth,2)
vb_dat$MG_rpkm <- outliers(vb_dat$MG_rpkm,2)

####################################################################
## ASSEMBLY STATISTICS AND VISUALIZATIONS
####################################################################
print("Begin visualizing data")
## Output filtering statistics table (text and tsv)

## Output assembly statistics table (text and tsv)
print("Print assembly statistics table")
sink(name_plot("assembly_stats.html"))
print(xtable(assembly.stats, html.table.attributes=""), type = "html")
sink()

write.table(assembly.stats, name_plot("assembly_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

## Output mapping stats table
print("Print metagenomic mapping statistics table")
sink(name_plot("mg_mapping_stats.html"))
print(xtable(MG.map.summary, html.table.attributes=""), type="html")
sink()

write.table(MG.map.summary, name_plot("mg_mapping_stats.txt"),
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

print("OVER")
####################################################################
## MAPPING STATISTICS AND VISUALIZATIONS
####################################################################

## Plot mappable reads density
print("Generating mapped reads plot")
var1 <-log10(c(all.dat$MG_reads,all.dat$MG_rpkm))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MG",nrow(all.dat)))
MG_mapped_reads<-data.frame(var1,var2)

png(name_plot("IMP-mg_reads_density.png"), width=350, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MG_mapped_reads,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Mappable reads", ylab=expression(log[10]*~"count"))
legend("bottomleft", fill =c("blue", "red"), legend = c("No of reads mapped", "RPKM normalized"))
dev.off()

## Plot coverage density
print("Generating MG coverage plot")
var1 <-c(all.dat$MG_cov,all.dat$MG_depth)
var2 <- c(rep("MG",nrow(all.dat)),rep("MG",nrow(all.dat)))
MG_coverage<-data.frame(var1,var2)

png(name_plot("IMP-mg_coverage_density.png"), width=350, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MG_coverage,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Coverage", ylab="fraction")
legend("bottomleft", fill =c("blue", "red"), legend = c("Coverage", "depth"))
dev.off()

## Plot vizbin scatter plot with length and MG coverage info
print("Generating vizbin plot for metagenomic coverage")
png(name_plot("IMP-mg_vizbin_length_cov.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MG_cov, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "coverage"))))
      ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MG depth info
print("Generating vizbin plot for metagenomic depth")
png(name_plot("IMP-mg_vizbin_length_depth.png"),width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MG_depth, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "depth"))))
      ) +
theme_nothing()
dev.off()

####################################################################
## VARIANT STATISTICS AND VISUALIZATIONS
####################################################################
## Plot variant count
var1 <-log10(c(all.dat$MG_var,all.dat$MG_dens))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MG",nrow(all.dat)))
MG_variant_count<-data.frame(var1,var2)

print("Generating variant count plots")
png(name_plot("IMP-mg_var_count.png") ,width=350, height=700)

par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MG_variant_count,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="MG variant (SNPs & INDELS)", ylab=expression(log[10]*~count))
legend("bottomleft", fill =c("blue", "red"), legend = c("No. of variants", "Variant density"))
dev.off()

## Plot vizbin scatter plot with length and MG variant info
# Create label
MG_var_label <- expression(bold(frac(variants[MG]/kb, MG[rpkm])))

print("Generating vizbin plot for metagenomic variant density")
png(name_plot("IMP-mg_vizbin_length_vardens.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(colour=MG_var_dens, size=log10(length), order=MG_var_dens), alpha=0.75) +
scale_colour_gradient(high="black", low="cyan") +
       guides(size=guide_legend(title=log10len),
	      colour=guide_colourbar(title=MG_var_label)
       ) +
theme_nothing()
dev.off()

## Save the R workspace
save.image(name_plot("mg_results.Rdat"))
