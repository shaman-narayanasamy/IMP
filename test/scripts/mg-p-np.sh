#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e

echo 'MG skip preprocessing preprocessing'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-np --current "snakemake preprocessing.done"
exit 0
