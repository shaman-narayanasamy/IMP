# Base container with tools needed for the IMP pipeline
#
# VERSION       1.2.0

FROM docker-r3lab.uni.lu/imp/imp-deps:1.2
MAINTAINER yohan.jarosz@uni.lu


######################################
#        Tools installations         #
######################################


## fastuniq
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/FastUniq-1.1.tar.gz -nv\
    && tar -xzf FastUniq-1.1.tar.gz \
    && cd FastUniq/source \
    && make \
    && mv fastuniq /usr/bin/. \
    && cd ../.. && rm -rf FastUniq* \
## Trimmomatic
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/Trimmomatic-Src-0.32.zip -nv\
    && unzip Trimmomatic-Src-0.32.zip \
    && cd trimmomatic-0.32 \
    && ant \
    && cp dist/jar/trimmomatic-0.32.jar /home/imp/lib/. \
    && cd .. && rm -rf *rimmomatic*

## idba ud
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/idba-1.1.1.tar.gz -nv\
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
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/cap3.linux.x86_64.tar -nv\
    && tar -xf cap3.linux.x86_64.tar \
    && cp CAP3/cap3 /usr/bin \
    && rm -rf cap3* CAP3 \
## bwa
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/bwa-0.7.9a.tar.bz2 -nv\
    && tar -xjf bwa-0.7.9a.tar.bz2 \
    && cd bwa-0.7.9a \
    && make \
    && mv bwa /usr/bin \
    && cd .. && rm -rf bwa*

## htsjdk - dependency of picard tools and FastQC (below)
RUN cd /home/imp/lib \
    && mkdir fastqc \
    && cd fastqc \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/htsjdk-1.138.tgz -O htsjdk-1.138.tgz -nv \
    && tar -xzf htsjdk-1.138.tgz \
    && cd htsjdk-1.138 \
    && ant htsjdk-jar \
    && cd .. && rm  htsjdk-1.138.tgz

## Picard tools - dependency of fastqc (below)
RUN cd /home/imp/lib \
    && cd /home/imp/lib/fastqc \
    &&  wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/picard.1.138.tgz -O picard-1.138.tgz -nv \
    && tar -xzf picard-1.138.tgz \
    && mv /home/imp/lib/fastqc/htsjdk-1.138  /home/imp/lib/fastqc/picard-1.138/htsjdk \
    && cd picard-1.138/ \
    && ant -lib lib/ant package-commands \
    && cd .. && rm  picard-1.138.tgz

## FastQC
RUN cd /home/imp/lib \
    && cd /home/imp/lib/fastqc \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/fastqc_v0.11.3.zip -nv \
    && unzip fastqc_v0.11.3.zip \
    && cd FastQC \
    && chmod 775 fastqc

## freebayes
RUN cd /home/imp/tmp \
    && git clone --recursive https://git-r3lab.uni.lu/R3/freebayes.git \
    && cd freebayes \
    && git checkout -b v0.9.16 \
    && make \
    && make install \
    && cd .. && rm -rf freebayes \
## vcftools
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/vcftools_0.1.12b.tar.gz -nv \
    && tar -xzf vcftools_0.1.12b.tar.gz \
    && cd vcftools_0.1.12b \
    && make \
    && make install \
    && cp -r bin/* /usr/bin \
    && cp -r perl/* /etc/perl/. \
    && cd .. && rm -rf vcftools*

## prokka
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/prokka-1.11.tar.gz -nv \
    && tar -xzf prokka-1.11.tar.gz \
    && cd prokka-1.11 \
    && cp bin/prokka* /usr/bin \
    && cp binaries/linux/* /usr/bin \
    && cd .. && rm -rf prokka* \
## parallel
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/parallel-20140422.tar.bz2 -nv \
    && tar -xjf parallel-20140422.tar.bz2 \
    && cd parallel-20140422 \
    && ./configure \
    && make \
    && make install \
    && cd .. && rm -rf parallel*

## sortmerna
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/sortmerna.2.0.tgz -nv \
    && tar -xzf sortmerna.2.0.tgz \
    && cd sortmerna-2.0 \
    && sh build.sh \
    && mv sortmerna indexdb_rna /usr/bin/. \
    && mv scripts/*.sh /home/imp/lib/. \
    && cd .. && rm -rf sortmerna*
## bedtools2
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/bedtools-2.24.0.tar.gz -O bedtools-2.24.0.tar.gz -nv \
    && tar -xzf bedtools-2.24.0.tar.gz \
    && cd bedtools2 \
    && make


## Install KronaPlot
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/KronaTools-2.5.tar -nv \
    && tar -xvf KronaTools-2.5.tar \
    && cd KronaTools-2.5 \
    && perl install.pl


## htslib
RUN cd /home/imp/tmp \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/htslib-1.2.1.tar.bz2 -nv \
    && tar -jxvf htslib-1.2.1.tar.bz2 \
    && cd htslib-1.2.1 \
    && ./configure && make && make install \
    && cd .. && rm -rf htslib-1.2.1

## Platypus
RUN cd /home/imp/lib \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/Platypus-latest.tgz -nv \
    && tar -xvzf Platypus-latest.tgz \
    && cd Platypus_0.8.1 \
    && bash buildPlatypus.sh


RUN cd /home/imp/lib \
    && rm -rf megahit \
    && git clone https://github.com/voutcn/megahit.git \
    && cd megahit \
    && make \
    && mv megahit* /usr/bin/.

## Vizbin JAR
RUN cd /home/imp/lib \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/VizBin-dist.jar -O VizBin-dist.jar -nv

RUN cd /home/imp/lib

## Quast/metaQuast
RUN cd /home/imp/lib \
    && git config --global http.postBuffer 2M \
    && git config --global http.maxRequestBuffer 100M \
    && git config --global core.compression 0 \
    && wget --no-check-certificate https://webdav-r3lab.uni.lu/public/R3lab/IMP/quast-3.1.zip -O quast.zip -nv\
    && unzip quast.zip \
    && cd quast-release_3.1 \
    && ln -s $PWD/metaquast.py /usr/bin/metaquast

######################
# set up environment #
######################

## add snakemake completion
RUN echo "`snakemake --bash-completion`" >> ~/.bashrc \
## add LD_LIBRARY_PATH
    && echo "export LD_LIBRARY_PATH=/usr/local/lib" >> ~/.bashrc

#####################
# Ship IMP codebase #
#####################

RUN mkdir /code \
    && cd /code \
    && git clone --branch 1.2.0 --single-branch https://git-r3lab.uni.lu/shaman.narayanasamy/IMP.git


######################
# runtime parameters #
######################
ENV PATH /home/imp/lib/trinityrnaseq_r20140717:/home/imp/lib/Platypus_0.8.1:/home/imp/lib/megahit:/home/imp/lib/fastqc/FastQC:$PATH
ENV LD_LIBRARY_PATH /usr/local/lib/:/home/imp/lib/Platypus_0.8.1
VOLUME ["/data", "/output", "/databases", "/code"]
WORKDIR /code
CMD ["snakemake", "ALL"]
