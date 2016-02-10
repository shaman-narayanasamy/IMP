# IMP

The Integrated Meta-omic Pipeline (IMP) is developed to perform large-scale, reproducible and automated integrative analysis of metagenomic and metatranscriptomic data. IMP also performs single omic (i.e. metagenomic-only and metatrancriptomic-only) analysis as an additional functionality.

IMP use various tools in order to perform an analysis (see `docs/DEPENDENCIES.md`). IMP does not provide an installation script for the tools it uses. It would be impossible to maintain.

Instead we provide the IMP workflow as is via the `snakemake` wrapper. Additionally, we provide a `docker container` with all the tools already installed inside.

Depending on your setup you might choose to use the workflow directly via `snakemake` or use the docker container via the `wrapper script`



# Documentation and website

All documentation and resources can be found : [here](http://r3lab.uni.lu/web/imp/doc.html)

All components used to develop IMP workflow care addressed under the [R3lab frozen pages](http://r3lab.uni.lu/frozen/imp/)
