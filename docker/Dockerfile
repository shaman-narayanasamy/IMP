# Base container with tools needed for the IMP pipeline
#
# VERSION       1.4.1

FROM docker-r3lab.uni.lu/imp/imp-deps:1.4.1
MAINTAINER yohan.jarosz@uni.lu

USER root


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


RUN mkdir -p /home/imp/code /home/imp/output /home/imp/data /home/imp/bin /home/imp/databases \
    && cd /home/imp/code \
    && git clone --branch 1.3 --single-branch https://git-r3lab.uni.lu/IMP/IMP.git . \
    && cd /usr && ln -s /home/imp/databases db


# R with checkpoint libraries
ADD dependencies.R /home/imp/lib/
RUN mkdir /root/.checkpoint \
    && echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > /root/.Rprofile \
    && Rscript -e "install.packages('checkpoint')" \
    && Rscript -e "library(checkpoint);checkpoint('2015-04-27', project='/home/imp/lib')" \
    && Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('genomeIntervals')"\
    && chmod -R 777 /root/.checkpoint /root/.Rprofile

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
