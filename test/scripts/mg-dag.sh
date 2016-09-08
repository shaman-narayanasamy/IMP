#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MG default DAGs'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg --current 'bash -c "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg --current 'bash -c "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg --current 'bash -c "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg --current 'bash -c "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report.pdf"'

echo 'MG no filtering DAGs'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-nf --current 'bash -c "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing-nf.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-nf --current 'bash -c "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly-nf.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-nf --current 'bash -c "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning-nf.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-nf --current 'bash -c "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report-nf.pdf"'

echo 'MG skip preprocessing DAGs'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-np --current 'bash -c "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing-np.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-np --current 'bash -c "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly-np.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-np --current 'bash -c "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning-np.pdf"'
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mg-np --current 'bash -c "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report-np.pdf"'
exit 0
