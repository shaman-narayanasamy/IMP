#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MG default DAGs'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg --current "snakemake -n"

echo 'MG no filtering DAGs'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-nf --current "snakemake -n"

echo 'MG skip preprocessing DAGs'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-np --current "snakemake -n"
exit 0
