# IMP

The Integrated Meta-omic Pipeline (IMP) is developed to perform large-scale, reproducible and automated integrative analysis of metagenomic and metatranscriptomic data. IMP also performs single omic (i.e. metagenomic-only and metatrancriptomic-only) analysis as an additional functionality.

IMP use various tools in order to perform an analysis (see `docs/DEPENDENCIES.md`). IMP does not provide an installation script for the tools it uses. It would be impossible to maintain.

Instead we provide the IMP workflow as is via the `snakemake` wrapper. Additionally, we provide a `docker container` with all the tools already installed inside.

Depending on your setup you might choose to use the workflow directly via `snakemake` or use the docker container via the `wrapper script`

* **wrapper script**: Everything is self-contained in a single docker container. Only 3 dependencies are needed in order to run IMP.
> Best on personal computers or when working in a mutually trusted environment.

* **snakemake**: Only the workflow is provided. You need to install correctly every tool on your machine.
> Best for cluster environment or when you don't have administrative rights to run docker.

Each run mode is respectively explained in the next sections.

---

### Wrapper script

In the same directory of this README file you can find `IMP` script.
This wrapper script will parse the command line arguments you provide
and translate everything to the docker container.


#### Dependencies

* [docker](https://docs.docker.com/installation)
* [git](http://www.git-scm.com)
* [python3]( https://www.python.org/downloads)

#### Usage

First you need to run

```bash
./IMP --init
```

to initialize and download databases needed to run IMP. It is only needed once. Then you could

```bash
./IMP -h  
```
to get some help.

---

### Snakemake

#### Dependencies

See `docs/DEPENDENCIES.md`.
You could look at the `docker/Dockerfile` to look how to install some tools.
We do not provide an installation script as it my vary a lot depending on the environment you are (cluster usage, module load, ...)


#### Usage

First you need to run

```bash
snakemake -s rules/init
```

to initialize and download databases needed to run IMP. It is only needed once. Then you could

```bash
snakemake -l # Get a list of available steps.
snakemake    # Launch the analysis
```

#### Configuration

We use a config file to pass variables to snakemake wrapper script. The default parameters are visible in `src/config.imp.json`.  You could override parameter via the file `conf/userconfig.imp.json`.

> Please do not override parameters directly on `src/config.imp.json` as it may be overridden with the next IMP update.

Eventually you could pass a different location for the config file:

```bash
CONFIGFILE=/home/imp/myconfigfile.json snakemake
```

Some parameters that can vary a lot can also be overridden via the command line
and take precedence over the config file:
* OUTPUTDIR
* MG
* MT
* DBPATH
* TMPDIR
* MEMTOTAL
* MEMCORE

You can provide them like this:

```bash
OUTPUTDIR=/home/imp/output MG="/home/imp/input/MG.R1.fq MG=/home/imp/input/MG.R2.fq" TMPDIR=/tmp snakemake
```

To see the full list of parameters and a description, see `docs/PARAMETERS.md`.
