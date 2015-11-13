### Quick Install

* Download (or clone with git) the latest IMP version:
[Versions available](https://git-r3lab.uni.lu/shaman.narayanasamy/IMP/tags)

* Install some dependencies:
    * [docker](https://docs.docker.com/installation)
    * [python3]( https://www.python.org/downloads)


* Untar/Unzip the folder and go inside the directory with your terminal.

* Start using it, you are ready ;)

~~~
./IMP --init  # initialize the workflow
./IMP -h      # get some help
~~~

### Full install

IMP uses the [Snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) workflow manager to handle
the analysis steps properly. The full installation is appropriate when you don't want to use the docker container
or you cannot because you are not in a mutually trusted environment
(see [Docker security documentation](https://docs.docker.com/articles/security/)) for more information about this.


As the installation step may vary a lot from one environment to another, you will have to install every tools manually
You could look at the **docker files** in the source code to look how to install some of the tools
(under `docker` directory.)



#### List of dependencies


* [openjdk](http://openjdk.java.net/) - 7
* [bioperl](http://www.bioperl.org/wiki/Main_Page) - 1.6.923-1
* [tabix](http://www.htslib.org/doc/tabix.html) - 0.2.6-2
* [fastuniq](http://sourceforge.net/projects/fastuniq/) - 1.1
* [samtools](http://samtools.sourceforge.net/)
* [gnuplot](http://www.gnuplot.info/)
* [python](http://python.org/) - 3.4
* [R](https://www.r-project.org/)
* [python](http://python.org/) - 2.7
* [gfortran](http://gcc.gnu.org/fortran/)
* [libatlas](http://math-atlas.sourceforge.net/)
* [docopt](http://docopt.org/)
* [bioservices](https://pypi.python.org/pypi/bioservices) - 1.3.5
* [numpy](http://www.numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [sklearn](http://scikit-learn.org/stable/index.html)
* [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) - 0.32
* [idba ud](http://wiki.hpc.ufl.edu/doc/IDBA-UD) - 1.1.1
* [cap3](http://seq.cs.iastate.edu/cap3.html)
* [bwa](http://bio-bwa.sourceforge.net/) - 0.7.9a
* [htsjdk](https://github.com/samtools/htsjdk) - 1.138
* [Picard tools](https://github.com/broadinstitute/picard)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - v0.11.3
* [Freebayes](https://github.com/ekg/freebayes) - v0.9.16
* [vcftools](http://vcftools.sourceforge.net/) - 0.1.12b
* [prokka](https://github.com/tseemann/prokka) - 1.11
* [parallel](https://www.gnu.org/software/parallel/) - 20140422
* [sortmerna](http://bioinfo.lifl.fr/RNA/sortmerna/) - 2.0
* [bedtools2](https://github.com/arq5x/bedtools2) - 2.24.0
* [KronaTools](https://github.com/marbl/Krona/wiki) - 2.5
* [htslib](http://www.htslib.org/) - 1.2.1
* [Platypus](https://github.com/andyrimmer/Platypus) - 0.8.1
* [megahit](https://github.com/voutcn/megahit)
* [Vizbin](https://github.com/claczny/VizBin)
* [Quast/metaQuast](http://bioinf.spbau.ru/en/metaquast) - 3.1

<div class="alert alert-info" role="alert">
    Some lower dependencies are not listed.
    For a complete list of all dependencies, please see what is installed in `docker/Dockerfile-dependencies`
</div>

### R dependencies

* genomeIntervals from [biocLite](http://bioconductor.org/install/).

We use the [checkpoint](http://projects.revolutionanalytics.com/documents/rrt/rrtpkgs/) library set to the `2015-04-27`
to install the following R packages:

* ggplot2
* gtools
* data.table
* reshape
* grid
* grDevices
* genomeIntervals
* stringr
* xtable
* beanplot
* psych
