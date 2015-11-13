# Introduction

IMP uses the [Snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) workflow manager to handle
the analysis steps properly.

You can execute the workflow directly on your computer/cluster using the `snakemake` command but you will have to install
each tools properly by yourself (see [Installation](installation.html)).

Alternatively, we also provide the workflow with every tools already installed within a `docker container`.
All you have to do is use the **ÌMP wrapper script**

Depending on your setup you might choose one or the other:

* **wrapper script**: Everything is self-contained in a single docker container. Only 3 dependencies are needed in order to run IMP.

> Best on personal computers or when working in a mutually trusted environment.

* **snakemake**: Only the workflow is provided. You need to install correctly every tool on your machine.

> Best for cluster environment or when you don't have administrative rights to run docker.


# Usage

IMP takes Meta-omic data as input. Understand paired fastq files.

In the following, Paired Metagenomic files are MG.r1.fq and MG.r2.fq .
Metatranscriptomic files are MT.r1.fq and MT.r2.fq.


## Wrapper script

In the root directory of IMP project you can find `IMP` script.
This wrapper script will parse the command line arguments you provide
and translate everything to the underlying **docker container**.

Each time a IMP command is run, it prints out the docker command used. You could also use the latter one if you
understand how docker works. It could provides you more control but **IMP** command line should be sufficient.


~~~
 # get help
./IMP -h
# Initialisation of the databases, on the first run only
./IMP --init    
# You could specify an another path
IMP --init -d /path/to/databases
# simple run with default options
./IMP -m input/mg.r1.fq -m input/mg.r2.fq -t input/mt.r1.fq -t input/mt.r2.fq -o output_directory
# alternative path for databases set in init mode
./IMP -m input/mg.r1.fq -m input/mg.r2.fq -t input/mt.r1.fq -t input/mt.r2.fq -o output_directory -d /path/to/databases
~~~

## Snakemake

You could use the workflow directly if you have all tools already installed (see [Installation](INSTALL.md)) or **enter** the container to use Snakemake directly:

~~~
./IMP -m input/mg.r1.fq -m input/mg.r2.fq -t input/mt.r1.fq -t input/mt.r2.fq -o output_directory --enter
~~~


~~~
# list all steps
snakemake -l
# initialise
snakemake -s rules/init
# launch the analysis
snakemake
~~~


# Configuration

IMP uses a config file to pass variables to the underlying Snakemake wrapper script.
See [Configuration](Configuration.md)
