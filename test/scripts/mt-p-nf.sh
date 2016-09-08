#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e

echo 'MT no filtering preprocessing'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-nf --current "snakemake preprocessing.done"

exit 0
