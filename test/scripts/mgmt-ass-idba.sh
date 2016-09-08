#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MGMT Assembly with IDBA default DAGs'
cp /mnt/data/input/*.preprocessed.fq $CI_PROJECT_DIR/output-mgmt/.
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -a idba -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mgmt --current "snakemake assembly.done"

exit 0
