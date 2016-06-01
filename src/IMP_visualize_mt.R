#!/bin/R

###################################################################################################
## Load required packages
###################################################################################################
print("Loading required R libraries")
require(genomeIntervals)

require(checkpoint)
checkpoint('2015-04-27', scanForPackages=FALSE, checkpointLocation="/root/")

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
MT.read.count_file  <- args[2]
MT.map.summary_file <- args[3]
MT.cov_file	    <- args[4]
MT.depth_file	    <- args[5]
MT.var_file	    <- args[6]
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

## Read counts MT
print("Read in MT read count file")
MT.read.count <- read.table(MT.read.count_file)
MT.read.count <- as.data.frame(cbind(as.character(MT.read.count$V2)[-nrow(MT.read.count)], MT.read.count$V1[-nrow(MT.read.count)]))
colnames(MT.read.count) <- c("file", "count")
MT.read.count$file <- as.character(MT.read.count$file)
# Divide by four for total sequences because original file is a line count
MT.read.count$count <- as.numeric(as.character(MT.read.count$count))/4

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

## MT coverage (samtools idxstat output)
print("Read in MT coverage file")
MT.cov <- read.table(MT.cov_file, colClasses=c("factor", rep("numeric", 6 )))[,c(1,4:7)]
colnames(MT.cov) <- c("contig", "MT_reads", "MT_bases_covered", "length", "MT_cov")

## MT depth (bedtools coverageBed output)
print("Read in MT depth file")
MT.depth <- read.table(MT.depth_file, colClasses=c("factor", "numeric"))
colnames(MT.depth) <- c("contig", "MT_depth")

## MT variant calls
print("Read in MT variant calls")
MT.var <- as.data.frame(table(read.table(MT.var_file)[,1]))
colnames(MT.var) <- c("contig", "MT_var")

## GC content file
print("Read in GC percentage file")
GC.dat <- read.delim(GC.dat_file)
colnames(GC.dat) <- c("contig", "GC")

### gff annotation file
print("Read in gff3 annotation file")
annot <- readZeroLengthFeaturesGff3(annot_file)

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

annot.4 <- merge(annot.2, annot.3, by="contig")

# vizbin points, the contig names are usually not included in this file,
# concatenate it before reading in
print("Reading in vizbin coordinates")
coords <- read.delim(coords_file, colClasses=c("factor", "numeric", "numeric"),
		     sep="\t", col.names=c("contig", "x", "y"), header=F)

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
all.dat <- merge(MT.cov, GC.dat, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MT.depth, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MT.var, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, annot.4, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, nucmer_res, by=c("contig"), all=T, incomparables=NA)
print("DONE: Merging data")

###################################################################################################
## Perform calculations and append it to the full table
###################################################################################################
# Calculate and merge data
print("Perform calculations")
save.image(name_plot("results.Rdat"))
# Get new column numbers
newcols <- ncol(all.dat) + 1
all.dat <- cbind(all.dat,
      gene_density(all.dat$no_of_genes, all.dat$length),
      coding_density(all.dat$total_gene_length, all.dat$length),
      contig_rpkm(all.dat$MT_reads, all.dat$length),
      var_density(all.dat$MT_var, all.dat$length, all.dat$MT_reads)
)
colnames(all.dat)[newcols:ncol(all.dat)] <- c(
			      "gene_dens",
			      "coding_dens",
			      "MT_rpkm",
			      "MT_var_dens"
			      )
# Get new column numbers
newcols <- ncol(all.dat) + 1

###################################################################################################
## Organize filtering statistics and create table
###################################################################################################
## Organize MT filtering procedures
print("Printing metagenomic filtering statistics")
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

# Remove outliers and infinite values
vb_dat$MT_var_dens <- outliers(vb_dat$MT_var_dens,2)
vb_dat$MT_depth <- outliers(vb_dat$MT_depth,2)
vb_dat$MT_rpkm <- outliers(vb_dat$MT_rpkm,2)

#write.table(vb_dat, "final.contig.merged_min1000_info_processed.txt", sep="\t", quote=F, row.names=F)

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

#print(assembly.stats_plot)

## Output mapping stats table
print("Print metagenomic mapping statistics table")
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

print("OVER")
####################################################################
## MAPPING STATISTICS AND VISUALIZATIONS
####################################################################

## Plot mappable reads density
print("Generating mapped reads plot")
var1 <-log10(c(all.dat$MT_reads,all.dat$MT_rpkm))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MT",nrow(all.dat)),rep("MT",nrow(all.dat)))
MT_mapped_reads<-data.frame(var1,var2) 

png(name_plot("IMP-MT_reads_density.png"), width=350, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MT_mapped_reads,  side = "both",log="auto", 
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Mappable reads", ylab=expression(log[10]*~"count"))
legend("bottomleft", fill =c("blue", "red"), legend = c("No of reads mapped", "RPKM normalized"))
dev.off()

## Plot coverage density
print("Generating MT coverage plot")
var1 <-c(all.dat$MT_cov,all.dat$MT_depth)
var2 <- c(rep("MT",nrow(all.dat)),rep("MT",nrow(all.dat)))
MT_coverage<-data.frame(var1,var2) 

png(name_plot("IMP-MT_coverage_density.png"), width=350, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MT_coverage,  side = "both",log="auto", 
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Coverage", ylab="fraction")
legend("bottomleft", fill =c("blue", "red"), legend = c("Coverage", "depth"))
dev.off()

## Plot vizbin scatter plot with length and MT coverage info
print("Generating vizbin plot for metagenomic coverage")
png(name_plot("IMP-MT_vizbin_length_cov.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MT_cov, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "coverage"))))
      ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MT depth info
print("Generating vizbin plot for metagenomic depth")
png(name_plot("IMP-MT_vizbin_length_depth.png"),width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MT_depth, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "depth"))))
      ) +
theme_nothing()
dev.off()

####################################################################
## VARIANT STATISTICS AND VISUALIZATIONS
####################################################################
## Plot variant count
var1 <-log10(c(all.dat$MT_var,all.dat$MT_dens))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MT",nrow(all.dat)),rep("MT",nrow(all.dat)))
MT_variant_count<-data.frame(var1,var2) 

print("Generating variant count plots")
png(name_plot("IMP-MT_var_count.png") ,width=350, height=700)

par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MT_variant_count,  side = "both",log="auto", 
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
main="MT variant (SNPs & INDELS)", ylab=expression(log[10]*~count))
legend("bottomleft", fill =c("blue", "red"), legend = c("No. of variants", "Variant density"))
dev.off()

## Plot vizbin scatter plot with length and MT variant info
# Create label
MT_var_label <- expression(bold(frac(variants[MT]/kb, MT[rpkm])))

print("Generating vizbin plot for metagenomic variant density")
png(name_plot("IMP-MT_vizbin_length_vardens.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(colour=MT_var_dens, size=log10(length), order=MT_var_dens), alpha=0.75) +
scale_colour_gradient(high="black", low="cyan") +
       guides(size=guide_legend(title=log10len),
	      colour=guide_colourbar(title=MT_var_label)
       ) +
theme_nothing()
dev.off()

## Save the R workspace
save.image(name_plot("results.Rdat"))


