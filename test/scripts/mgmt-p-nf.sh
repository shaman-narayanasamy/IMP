#!/usr/bin/env bash
cd $CI_PROJECT_DIR

set -e
echo 'MGMT no filtering preprocessing'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-nf --current "snakemake preprocessing.done"
exit 0
