#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MT default preprocessing'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt --current "snakemake preprocessing.done"

exit 0
