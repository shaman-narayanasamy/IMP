# How to build a new tarball version from docker files

> This README assumes that you are in the `docker` directory of the IMP project.

## Increment Version number

Edit the docker file `docker/Dockerfile`:

Increment the version number in that file.

> The version number can be anything: e.g '1.1.2' or 'my-feature'

Edit the `IMP` script and put the `LATEST_TAG` variable to the same version.

## Build the docker images locally

> In order to build the docker images, you must have [Docker](https://docs.docker.com/installation/) installed.

### Build dependencies

    docker build -t docker-r3lab.uni.lu/imp/imp-deps:1.1 -f Dockerfile-dependencies .

> 'docker-r3lab.uni.lu/imp/imp-deps:1.1' is the image name that we will have to give to the main Docker file.

### Build the main Docker file.

Edit the docker file and change `FROM` to the image name you gave in the previous step.

Then build the file

    docker build -t docker-r3lab.uni.lu/imp/imp:<version> .


## Build the tarball

    docker save docker-r3lab.uni.lu/imp/imp:<version> > imp-<version>.tar
    gzip -9 imp-<version>.tar


## Put the tarball online

    Must be under `https://webdav-r3lab.uni.lu/public/R3lab/IMP/dist/`