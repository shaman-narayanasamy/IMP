# Base container with tools needed for the IMP pipeline
#
# VERSION       1.4.1

FROM docker-r3lab.uni.lu/imp/imp-deps:1.4.1
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

## CheckM
RUN cd /home/imp/lib \
    && pip2.7 install pysam dendropy ScreamingBackpack \
    && wget https://github.com/Ecogenomics/CheckM/archive/v1.0.6.tar.gz \
    && tar xf v1.0.6.tar.gz \
    && cd CheckM-1.0.6/ \
    && python2.7 setup.py install \
    && cd .. \
    && rm -r CheckM-1.0.6/ \
    && wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Linux-v1.1.alpha17.zip \
    && unzip pplacer-Linux-v1.1.alpha17.zip \
    && cp pplacer-Linux-v1.1.alpha17/guppy /usr/bin \
    && cp pplacer-Linux-v1.1.alpha17/pplacer /usr/bin \
    && rm -r pplacer-Linux-v1.1.alpha17

##CheckM Data
RUN cd /home/imp/lib \
    && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_v1.0.6.tar.gz

RUN cd /home/imp/lib \
    && mkdir checkm_data_v1.0.6 \
    && tar -xf checkm_data_v1.0.6.tar.gz -C checkm_data_v1.0.6 \
    #this might open interactively:
    && echo checkm_data_v1.0.6 | checkm data setRoot checkm_data_v1.0.6/


##Amphora2
RUN cd /home/imp/lib \
    && git clone https://github.com/mherold1/AMPHORA2 \
    && wget https://github.com/stamatak/standard-RAxML/archive/v8.2.8.tar.gz -O stamatak-standard-RAxML.tar.gz \
    && tar xf stamatak-standard-RAxML.tar.gz \
    && cd standard-RAxML-8.2.8; make -f Makefile.gcc; make -f Makefile.PTHREADS.gcc; cp raxmlHPC* /usr/bin/. \
    && cd .. \
    && (echo y;echo o conf prerequisites_policy follow;echo o conf commit)|cpan -e install Parallel::ForkManager \
#    && wget ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz \
#    && tar -xf EMBOSS-6.5.7.tar.gz \
#    && cd EMBOSS-6.5.7 \
#    && ./configure \
#    && make \
#    && make install
    && apt-get install emboss

##pullseq
RUN cd /home/imp/lib \
   # && wget https://github.com/bcthomas/pullseq/releases/download/1.0.2/pullseq_v1.0.2_linux64.zip \
   # && unzip pullseq_v1.0.2_linux64.zip \
   # && rm seqdiff \
   # && mv pullseq /usr/bin
    && git clone https://github.com/bcthomas/pullseq \
    && cd pullseq \
    && ./bootstrap \
    && ./configure \
    && make \
    && make install

##Phylophlan
#RUN cd /home/imp/lib \
##    && wget https://bitbucket.org/nsegata/phylophlan/get/default.tar.gz \
##    && tar -xf default.tar.gz \
#    && pip install mercurial \
#    && hg clone https://bitbucket.org/nsegata/phylophlan \
##    && cd phylophlan \
##for development version: (stick to old one for now)
##    && hg pull && hg update dev \ 
#    && wget http://www.microbesonline.org/fasttree/FastTree \
#    && chmod +x FastTree \
#    && cp FastTree /usr/bin \
##usearch binary copied from gaia for now and put in imp code -> need to solve this 
#
#    #move usearch binary to webdav? this uses my downloadlink, from 2.5.2016 , not sure how long this works.
##    && wget --user-agent=Mozilla --content-disposition -E -c http://drive5.com/cgi-bin/upload3.py?license=2016050205235811143 \
# #   && chmod +x usearch8.1.1861_i86linux32 \
#  #  && cp usearch8.1.1861_i86linux32 /usr/bin/usearch \
##    && cd  \
##    && ln -s $PWD/phylophlan.py /usr/bin/phylophlan

######################
# set up environment #
######################

## add snakemake completion
RUN echo "`snakemake --bash-completion`" >> ~/.bashrc \
## add LD_LIBRARY_PATH
    && echo "export LD_LIBRARY_PATH=/usr/local/lib" >> ~/.bashrc


#### Add IMP user
RUN useradd -d /home/imp imp \
    && echo 'imp:imp' |chpasswd \
    && echo "imp ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/imp


#####################
# Ship IMP codebase #
#####################

USER root
RUN mkdir -p /home/imp/code /home/imp/output /home/imp/data /home/imp/bin /home/imp/databases \
    && cd /home/imp/code \
    && git clone --branch 1.3 --single-branch https://git-r3lab.uni.lu/IMP/IMP.git .


#### Add gosu
ENV GOSU_VERSION 1.7
RUN set -x \
    && apt-get update && apt-get install -y --no-install-recommends ca-certificates && rm -rf /var/lib/apt/lists/* \
    && wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$(dpkg --print-architecture)" \
    && wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$(dpkg --print-architecture).asc" \
    && export GNUPGHOME="$(mktemp -d)" \
    && gpg --keyserver ha.pool.sks-keyservers.net --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4 \
    && gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu \
    && rm -r "$GNUPGHOME" /usr/local/bin/gosu.asc \
    && chmod +x /usr/local/bin/gosu \
    && gosu nobody true

COPY entrypoint.sh /usr/local/bin/entrypoint.sh


RUN chown imp:imp -R /home/imp

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]



######################
# runtime parameters #
######################
ENV AMPHORA2_home /home/imp/lib/AMPHORA2  ###CHECK IF THIS IS CORRECT?
ENV PATH /home/imp/lib/pullseq/src:/home/imp/lib/trinityrnaseq_r20140717:/home/imp/lib/Platypus_0.8.1:/home/imp/lib/megahit:/home/imp/lib/fastqc/FastQC:$PATH
ENV LD_LIBRARY_PATH /usr/local/lib/:/home/imp/lib/Platypus_0.8.1
VOLUME ["/home/imp/data", "/home/imp/output", "/home/imp/databases", "/home/imp/code"]
WORKDIR /home/imp/code

CMD ["snakemake", "ALL"]
