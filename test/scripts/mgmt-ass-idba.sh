#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MGMT Assembly with IDBA default DAGs'
OUT=$CI_PROJECT_DIR/output-mgmt-idba
mkdir -p $OUT
cp /mnt/data/input/*.preprocessed.fq $OUT/.
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -a idba -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $OUT --current "snakemake assembly.done"

exit 0
