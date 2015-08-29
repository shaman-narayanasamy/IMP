#!/bin/R

###################################################################################################
## This file contains all the functions required for all the plots
## 
###################################################################################################

## calculate coverage using the "traditional" method
get_coverage=function(reads_mapped, length, read_len){
    M <- reads_mapped
    L <- length
    R <- read_len

    # coverage*contig length/read length
    C <- (M*L)/R
    return(C)
}

####################################################################
## Contig based "RPKM" values, similar to used in
## Muller et al. (2014, Nat. Comm.), ## only here it
## is applied to
## every contig
contig_rpkm=function(reads_mapped, length){
    N <- sum(reads_mapped)
    R <- reads_mapped
    L <- length

    # reads mapped/([length of contig]/1000)/([total reads]/10^6)
    R/(L/1000)/(N/10^6)
}

####################################################################
## Calculate variant density:
## i)  Traditional variats per kilo base
## ii) Normalized by contig rpkm (Muller et al., 2014)
var_density=function(variants, length, reads_mapped){
    V <- variants
    L <- length
    rpkmC <- contig_rpkm(reads_mapped, L)

    D <- (V/L)/1000/rpkmC
    return(D)
}

####################################################################
## Calculate gene density
gene_density=function(total_genes, length){
    G <- total_genes
    L <- length

    C <- (G)/(L/1000)
    return(C)
}
####################################################################
## Calculate coding density
coding_density=function(total_len_genes, length){
    G <- total_len_genes
    L <- length

    C <- (G/1000)/(L/1000)
    return(C)
}

####################################################################
## Calculate N50
get_n50=function(lengths){
    x=lengths
    x[cumsum(x) > sum(x)/2][1]
}


####################################################################
## Get assembly statistics
get_stats=function(dat){
    contigs <- nrow(dat)
    N50 <- get_n50(dat$length)
    max_len <- max(dat$length)
    mean_len <- mean(dat$length)
    med_len <- median(dat$length)
    #MG_mapped <- sum(dat$MG_mapped)
    #MT_mapped <- sum(dat$MT_mapped)
    total_length <- sum(dat$length)
    return(c(contigs,N50,max_len,mean_len,med_len,total_length))
}


####################################################################
## Function to scale data between 0 and 1 (for plotting)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

####################################################################
## Filter out outliers and set them as floor/ceiling value (min/max)
outliers=function(z, dist){
    z[is.infinite(z)] <- max(z[is.finite(z)])
    z[which(z > mean(z) + dist*sd(z))] = round(mean(z) + dist*sd(z))
    z[which(z < mean(z) - dist*sd(z))] = round(mean(z) - dist*sd(z))
    return(z)
}

####################################################################
## Function to name files
name_plot=function(name){
    filename <- paste(out_dir, name, sep='/')
    return(filename)
}

####################################################################
# Function to obtain number of character (char) occurences
# within string
countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}

####################################################################
# Function to output filename without the path
get_file_name <- function(file_path){
    n <- countCharOccurrences('/', as.character(file_path))
    filename <- str_split_fixed(file_path, "/", n+1)[n+1]
    return(filename)
}

####################################################################
# Function to obtain filtering information (contained in file names)
filtering <- function(file){
    n <- countCharOccurrences("\\.", as.character(file))
    filter <- str_split_fixed(file, "\\.", n+1)[n]
    return(filter)
}

## Function for white background theme with no axes
theme_nothing <- function(base_size = 12, base_family = "Helvetica")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
            rect             = element_blank(),
            line             = element_blank(),
            axis.ticks.margin = unit(0, "lines"),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
           )
}

## Function for black background theme with no axes
theme_black <- function(base_size = 12, base_family = "Helvetica")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
	    panel.background = element_rect(fill="black", colour="black"),
            line             = element_blank(),
            axis.ticks.margin = unit(0, "lines"),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
           )
}

## Function for gray background theme with no axes
theme_gray <- function(base_size = 12, base_family = "Helvetica")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
	    panel.background = element_rect(fill="gray25", colour="gray25"),
            line             = element_blank(),
            axis.ticks.margin = unit(0, "lines"),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
           )
}


## Set maximum value for plots based on standard deviation
set_max_sd=function(x, dist){
    max_val <- mean(x, na.rm=T) + dist*sd(x, na.rm=T)
    return(max_val)
}

## Set maximum value for plots based on percentage
set_max_perc=function(x, percentage){
    max_val <- max(x[is.finite(x)], na.rm=T)*(percentage/100)
}

## Special string for length
log10len <- expression(bold(atop("Contig length", paste("(",log[10], bp, ")"))))

## Metagenomic and metatranscriptomic labels
mgmt_labs <- c("metagenomic","metatranscriptomic")


