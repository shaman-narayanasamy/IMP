# Base container with tools needed for the IMP pipeline
#
# VERSION       0.1

FROM ubuntu:trusty
MAINTAINER yohan.jarosz@uni.lu


################
# dependencies #
################
RUN apt-get update && apt-get -yq upgrade \
    && apt-get install -yq libboost-dev \
    && apt-get install -yq make cmake wget build-essential unzip \
    && apt-get install -yq openjdk-7-jdk ant git tig vim \
    && apt-get install -yq libboost-filesystem1.54-dev libboost-program-options1.54-dev libboost-regex1.54-dev libboost-system1.54-dev libboost-thread1.54-dev libboost-system1.54-dev \
    && apt-get install -yq bioperl=1.6.923-1 tabix=0.2.6-2 zlib1g zlib1g-dev \
    && apt-get install -yq samtools gnuplot python3 \
    && apt-get install -yq r-base r-base-dev \
    && apt-get install -yq libncurses5 libncurses5-dev sudo libatlas-base-dev python2.7 gfortran python-dev \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --set java /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java


## Python dependencies and bioservices. Using version 1.3.5 because script doesn't work with latest version
RUN mkdir -p /home/imp/tmp /home/imp/lib \
    && cd /home/imp/tmp \
    && wget https://bootstrap.pypa.io/get-pip.py \
    && python3 get-pip.py \
    && python2.7 get-pip.py \
    && pip3 install snakemake docopt \
    && pip2.7 install numpy \
    && pip2.7 install scipy scikit-learn docopt \
    && pip install -Iv https://pypi.python.org/packages/source/b/bioservices/bioservices-1.3.5.tar.gz \
    && rm get-pip.py


## R
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile \
    && Rscript -e "install.packages('ggplot2')" \
    && Rscript -e "install.packages('gtools')" \
    && Rscript -e "install.packages('data.table')" \
    && Rscript -e "install.packages('reshape')" \
    && Rscript -e "install.packages('stringr')" \
    && Rscript -e "install.packages('xtable')" \
    && Rscript -e "install.packages('psych')" \
    && Rscript -e "install.packages('beanplot')" \
    && Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('genomeIntervals')"


######################################
#        Tools installations         #
######################################


## fastuniq
RUN cd /home/imp/tmp \
    && wget http://downloads.sourceforge.net/project/fastuniq/FastUniq-1.1.tar.gz \
    && tar -xzf FastUniq-1.1.tar.gz \
    && cd FastUniq/source \
    && make \
    && mv fastuniq /usr/bin/. \
    && cd ../.. && rm -rf FastUniq* \
## Trimmomatic
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.32.zip \
    && unzip Trimmomatic-Src-0.32.zip \
    && cd trimmomatic-0.32 \
    && ant \
    && cp dist/jar/trimmomatic-0.32.jar /home/imp/lib/. \
    #&& mkdir -p /usr/db/trimmomatic \
    #&& cp adapters/* /usr/db/trimmomatic/. \
    && cd .. && rm -rf *rimmomatic*

## idba ud
RUN cd /home/imp/tmp \
    && wget http://hku-idba.googlecode.com/files/idba-1.1.1.tar.gz \
    && tar -xzf idba-1.1.1.tar.gz \
    && cd idba-1.1.1 \
    && sed -i -e 's/static const uint32_t kMaxShortSequence = 128/static const uint32_t kMaxShortSequence = 2048/' src/sequence/short_sequence.h \
    && ./configure \
    && make \
    && make install \
    && mv bin/idba_ud /usr/bin \
    && mv bin/fq2fa /usr/bin \
    && cd .. && rm -rf idba* \
## cap3
    && wget http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar \
    && tar -xf cap3.linux.x86_64.tar \
    && cp CAP3/cap3 /usr/bin \
    && rm -rf cap3* CAP3 \
## bwa
    && wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.9a.tar.bz2 \
    && tar -xjf bwa-0.7.9a.tar.bz2 \
    && cd bwa-0.7.9a \
    && make \
    && mv bwa /usr/bin \
    && cd .. && rm -rf bwa*

## htqc toolbox
RUN cd /home/imp/tmp \
    && wget http://downloads.sourceforge.net/project/htqc/htqc-1.91.1-Source.tar.gz \
    && tar -xzf htqc-1.91.1-Source.tar.gz \
    && cd htqc-1.91.1-Source \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd ../.. && rm -rf htqc* \
    && cd .. && rm -rf cufflinks* \
## freebayes
    && git clone --recursive git://github.com/ekg/freebayes.git \
    && cd freebayes \
    && git checkout -b v0.9.16 v0.9.16 \
    && make \
    && make install \
    && cd .. && rm -rf freebayes \
## vcftools
    && wget http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz \
    && tar -xzf vcftools_0.1.12b.tar.gz \
    && cd vcftools_0.1.12b \
    && make \
    && make install \
    && cp -r bin/* /usr/bin \
    && cp -r perl/* /etc/perl/. \
    && cd .. && rm -rf vcftools*

## prokka
RUN cd /home/imp/tmp \
    && wget http://www.vicbioinformatics.com/prokka-1.11.tar.gz \
    && tar -xzf prokka-1.11.tar.gz \
    && cd prokka-1.11 \
    && cp bin/prokka* /usr/bin \
    && cp binaries/linux/* /usr/bin \
    && cd .. && rm -rf prokka* \
## parallel
    && wget ftp://ftp.gnu.org/gnu/parallel/parallel-20140422.tar.bz2 \
    && tar -xjf parallel-20140422.tar.bz2 \
    && cd parallel-20140422 \
    && ./configure \
    && make \
    && make install \
    && cd .. && rm -rf parallel* \
## sortmerna
    && wget https://github.com/biocore/sortmerna/archive/2.0.tar.gz -O sortmerna.2.0.tgz \
    && tar -xzf sortmerna.2.0.tgz \
    && cd sortmerna-2.0 \
    && sh build.sh \
    && mv sortmerna indexdb_rna /usr/bin/. \
    && mv scripts/*.sh /home/imp/lib/. \
    && cd .. && rm -rf sortmerna* \
## bedtools2
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 \
    && git checkout v2.22.0 \
    && make

## bh_tsne
RUN cd /home/imp/tmp \
    && git clone https://github.com/claczny/VizBin.git \
    && cd VizBin/src/backend/bh_tsne \
    && g++ quadtree.cpp tsne.cpp -o bh_tsne -O3 -I./CBLAS/include -L./ -lcblas \
    && mv bh_tsne /usr/bin/ \
    && cd /home/imp/tmp \
    && rm -rf VizBin

## Install KronaPlot
RUN cd /home/imp/tmp \
    && wget "http://downloads.sourceforge.net/project/krona/KronaTools (Mac, Linux)/KronaTools-2.5.tar" \
    && tar -xvf KronaTools-2.5.tar \
    && cd KronaTools-2.5 \
    && perl install.pl


## htslib
RUN cd /home/imp/tmp \
    && wget https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2 \
    && tar -jxvf htslib-1.2.1.tar.bz2 \
    && cd htslib-1.2.1 \
    && ./configure && make && make install \
    && cd .. && rm -rf htslib-1.2.1

## Platypus
RUN cd /home/imp/lib \
    && wget http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz \
    && tar -xvzf Platypus-latest.tgz \
    && cd Platypus_0.8.1 \
    && bash buildPlatypus.sh

## Megahit
RUN cd /home/imp/lib \
    && git clone https://github.com/voutcn/megahit.git \
    && cd megahit \
    && git checkout v0.2.1 \
    && make \
    && mv megahit* /usr/bin/. \
    && mv sdbg_builder_cpu /usr/bin/.


######################
# ship pipeline code #
######################
# (commented for now)
#RUN cd /home/imp \
#    && git clone https://git-r3lab.uni.lu/shaman.narayanasamy/IMP.git

######################
# set up environment #
######################

## add snakemake completion
RUN echo "`snakemake --bash-completion`" >> ~/.bashrc \
## add LD_LIBRARY_PATH
    && echo "export LD_LIBRARY_PATH=/usr/local/lib" >> ~/.bashrc

######################
# runtime parameters #
######################
ENV PATH /home/imp/lib/trinityrnaseq_r20140717:/home/imp/lib/Platypus_0.8.1:/home/imp/lib/megahit:$PATH

VOLUME ["/shared"]
WORKDIR /home/imp/integrated-metaomic-pipeline
CMD ["/bin/bash"]
