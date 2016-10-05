#!/bin/R

require(checkpoint)
checkpoint('2016-06-20', scanForPackages=FALSE, checkpointLocation="/home/imp/lib", project="/home/imp/lib")

require(stringr)
require(ggplot2)
require(RColorBrewer)

print("START: Reading arguments")
args <- commandArgs(trailingOnly = TRUE)
vb.points.file <- args[1]
mb.summary.file <- args[2]
mb.contigs.file <- args[3]
contig.len.file <- args[4]
function_script <- args[5]
print("DONE: Reading arguments")

## Load IMP visualization functions
print("START: Reading functions")
source(function_script)
print("DONE: Reading functions")

## Read in VizBin points
print("Read in VizBin points")
vb.points <- read.table(vb.points.file)
colnames(vb.points) <- c("contig", "x", "y")

## Read in MaxBin summary file
print("Read in MaxBin summary table")
mb.summary <- read.table(mb.summary.file, sep = "\t", header = T)
colnames(mb.summary) <- c("bin.name", "abundance", "completeness", "genome.size", "CG.content")
mb.summary$completeness <- as.numeric(gsub("%", "", mb.summary$completeness))

## Read in maxbin contigs2bin mappings
print("Read in MaxBin contigs to bin mappings")
mb.contigs <- read.table(mb.contigs.file)
colnames(mb.contigs) <- c("contig", "bin.name")

## Read in contig length file
print("Read in contig lenghts")
contig.len <- read.table(contig.len.file)
colnames(contig.len) <- c("contig", "length")

## Merge data
print("Merging MaxBin and VizBin data")
vb.mb.dat <- merge(vb.points, mb.contigs, all.x=T, incomparables="NA")
vb.mb.dat <- merge(vb.mb.dat, mb.summary, all.x=T, incomparables="NA")
vb.mb.dat <- merge(vb.mb.dat, contig.len, all.x=T, incomparables="NA")

print("Generate colour scheme")
essPal <- colorRampPalette(brewer.pal(11,"Spectral"))(111)[111:1]

mb.vb.plot <- ggplot(vb.mb.dat, aes(x=x, y=y)) + 
geom_point(aes(size=log10(length), alpha=abundance, colour=completeness)) +
scale_colour_gradientn(colours=essPal,
                       na.value="gray75",
                       guide="colourbar",
                       guide_legend(title="Completeness (%)")
                       ) +
scale_alpha_continuous(range = c(0.5, 1),
                       guide_legend(title="Abundance")
                       ) +
guides(size=guide_legend(title=log10len, order=1)) +
theme_nothing()

print("START: Visualization")
png("Binning/MaxBin/IMP-MaxBin-vizbin_length_completeness_abundance.png", width=plotWidth, height=plotHeight)
mb.vb.plot
dev.off()
print("DONE: Visualization")
