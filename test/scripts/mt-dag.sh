#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MT default DAGs'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt --current 'bash -c "snakemake --dag"'

echo 'MT no filtering DAGs'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-nf --current 'bash -c "snakemake --dag"'

echo 'MT skip preprocessing DAGs'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-np --current 'bash -c "snakemake --dag"'
exit 0
