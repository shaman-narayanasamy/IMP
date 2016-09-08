#!/usr/bin/env bash
cd $CI_PROJECT_DIR
set -e
echo 'MT default DAGs'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt --current 'bash -c "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt --current 'bash -c "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt --current 'bash -c "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt --current 'bash -c "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report.pdf"'

echo 'MT no filtering DAGs'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-nf --current 'bash -c "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing-nf.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-nf --current 'bash -c "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly-nf.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-nf --current 'bash -c "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning-nf.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-nf --current 'bash -c "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report-nf.pdf"'

echo 'MT skip preprocessing DAGs'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-np --current 'bash -c "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing-np.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-np --current 'bash -c "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly-np.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-np --current 'bash -c "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning-np.pdf"'
./IMP -t $CI_PROJECT_DIR/test/MT.R1.small.fq -t $CI_PROJECT_DIR/test/MT.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -d /mnt/data/db -o $CI_PROJECT_DIR/output-mt-np --current 'bash -c "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report-np.pdf"'
