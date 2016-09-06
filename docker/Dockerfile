# Base container with tools needed for the IMP pipeline
#
# VERSION       1.4

FROM docker-r3lab.uni.lu/imp/imp-tools:1.4
MAINTAINER yohan.jarosz@uni.lu

USER root


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
RUN mkdir /home/imp/lib/.checkpoint \
    && echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > /home/imp/.Rprofile \
    && Rscript -e "install.packages('checkpoint', repos='https://cloud.r-project.org/')" \
    && Rscript -e "library(checkpoint);checkpoint('2015-04-27', project='/home/imp/lib', checkpointLocation='/home/imp/lib')" \
    && Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('genomeIntervals')"\
    && chmod -R 777 /home/imp/lib/.checkpoint /home/imp/.Rprofile

# add Snakemake
RUN pip3 install snakemake==3.7.1

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
ENV PATH /home/imp/lib/FragGeneScan1.30:/home/imp/lib/MaxBin-2.2.1:/home/imp/lib/pullseq/src:/home/imp/lib/Platypus_0.8.1:/home/imp/lib/megahit:/home/imp/lib/fastqc/FastQC:$PATH
ENV LD_LIBRARY_PATH /usr/local/lib/:/home/imp/lib/Platypus_0.8.1
VOLUME ["/home/imp/data", "/home/imp/output", "/home/imp/databases", "/home/imp/code"]
WORKDIR /home/imp/code

CMD ["snakemake", "ALL"]
