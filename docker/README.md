# Docker How To

**IMP** is provided as a docker container with all dependencies shipped inside.
You 'only' need to install docker and you are ready to use **IMP**.

For more information about `docker` and how to install and use it, please refers to the [Official documentation](https://docs.docker.com/).


## Run IMP latest container

* Run latest version of **IMP**. For the first run, it will download all dependencies. It can take some time.
```bash
docker run imp:latest
... downloading data ...
Error: Snakefile "Snakefile" not present.
```

## Run IMP on your data

* The complete command is

```bash
docker run \
-v $DATA_DIR:/data \
-v $OUTPUT_DIR:/output \
-e MG="/data/$MGR1 /data/$MGR2" \
-e MT="/data/$MTR1 /data/$MTR2" \
imp:latest
```
* `docker run`: Run a docker container

* `-v $DATA_DIR:/data` : Where `$DATA_DIR` is the directory under your data is located. (*All paths must be absolute**)

>You need to provide the container your data. By default, **IMP** expect that you have metagenomic and metatransciptiomic data and will look into the container `/data` directory. You do it with the `-v` (for volume) flag.

* `-v OUTPUT_DIR:/output`: Where `$OUTPUT_DIR` is the directory where output data will be located.

> As the `$DATA_DIR`, you need to mount an output directory to hold all result files from **IMP**.

* `-e MG="/data/$MGR1 /data/$MGR2"`: Where `$MGR1` and `$MGR2` are the metagenomic paired-end reads.

* `-e MT="/data/$MTR1 /data/$MTR2"`: Where `$MTR1` and `$MTR2` are the metatransciptiomic paired-end reads.

* `imp:latest`: The name of the container to use (<name>:<version>)

### Example


```bash
$ pwd
/Users/imp/data
$ ls
MG.R1.fq MG.R2.fq MT.R1.fq MT.R2.fq
$ ls /Users/imp/build

$ docker run \
-v /Users/imp/data:/data \
-v /Users/imp/build:/output \
-e MG="/data/MG.R1.fq /data/MG.R2.fq" \
-e MT="/data/MT.R1.fq MG=/data/MT.R2.fq" \
imp:latest
$ ls /Users/imp/build
IMP.html
Preprocessing
Assembly
Mapping
Analysis
...
```

## Use a configuration file.
Many parameters can be overridden in `IMP` via a config file. You may provide using it with the `-e` flag and using the `CONFIGFILE` environment variable.

```bash
docker run \
-e CONFIGFILE="path/to/config/file.json"
 ... other parameters ...
```
## Build new container image

```bash
cd <IMP_DIR>/docker
docker build -t <name> .
```
Name can be abything you want. Then you can

```bash
docker run <name>
```

## Use development version of IMP

You need to *erase* imp source code in the container and mount the latest version of the source code.


```bash
$ cd $DIR
$ git clone https://git-r3lab.uni.lu/shaman.narayanasamy/IMP.git
$ cd IMP
$ git checkout dev
$ docker run \
-v $DIR/IMP:/home/imp/integrated-metaomic-pipeline \
 ... other parameters ...

```
Any change to the source code will be applied immediatly inside the container.

## Log in inside the container

If you want to enter the container and use `IMP` from there, you'll have to erase the docker entrypoint:

```bash
docker run \
--entrypoint /bin/bash -it \
 ... other parameters here ...
```
