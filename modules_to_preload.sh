#!/bin/bash -l 

#source this file before execution of snakefile
module load Python
source /mnt/nfs/projects/ecosystem_biology/local_tools/IMP/bin/activate

module load MEGAHIT
module load BWA 
module load SAMtools 
module load BEDTools
module load OpenBLAS 
module load Boost/1.53.0-ictce-5.3.0 

export PATH=$PATH:/mnt/nfs/projects/ecosystem_biology/local_tools/idba-1.1.1.icc/bin

module load CAP3

#symbolic links for prokka db 
module load prokka

export PATH=$PATH:/mnt/nfs/projects/ecosystem_biology/local_tools/tabix-0.2.6
export PATH=$PATH:/mnt/nfs/projects/ecosystem_biology/local_tools/gkno_launcher/tools/freebayes/bin
export PATH=$PATH:/mnt/nfs/projects/ecosystem_biology/local_tools/vcftools/bin
export PERL5LIB=$PERL5LIB:/mnt/nfs/projects/ecosystem_biology/local_tools/vcftools/perl
export PATH=$PATH:/mnt/nfs/projects/ecosystem_biology/local_tools/Platypus/Platypus_0.7.9.1

module load R
Rscript -e "install.packages('beanplot')"

module list
