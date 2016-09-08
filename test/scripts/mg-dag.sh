#!/usr/bin/env bash
cd $CI_PROJECT_DIR
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -o $CI_PROJECT_DIR/output-mg --current "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -o $CI_PROJECT_DIR/output-mg --current "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -o $CI_PROJECT_DIR/output-mg --current "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/default.conf.json -o $CI_PROJECT_DIR/output-mg --current "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report.pdf"

./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -o $CI_PROJECT_DIR/output-mg-nf --current "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing-nf.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -o $CI_PROJECT_DIR/output-mg-nf --current "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly-nf.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -o $CI_PROJECT_DIR/output-mg-nf --current "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning-nf.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_filtering.conf.json -o $CI_PROJECT_DIR/output-mg-nf --current "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report-nf.pdf"

./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -o $CI_PROJECT_DIR/output-mg-np --current "snakemake --dag preprocessing.done | dot -Tpdf > $CI_PROJECT_DIR/preprocessing-np.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -o $CI_PROJECT_DIR/output-mg-np --current "snakemake --dag assembly.done | dot -Tpdf > $CI_PROJECT_DIR/assembly-np.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -o $CI_PROJECT_DIR/output-mg-np --current "snakemake --dag binning.done | dot -Tpdf > $CI_PROJECT_DIR/binning-np.pdf"
./IMP -m $CI_PROJECT_DIR/test/MG.R1.small.fq -m $CI_PROJECT_DIR/test/MG.R2.small.fq -c $CI_PROJECT_DIR/test/no_preprocessing.conf.json -o $CI_PROJECT_DIR/output-mg-np --current "snakemake --dag report.done | dot -Tpdf > $CI_PROJECT_DIR/report-np.pdf"
