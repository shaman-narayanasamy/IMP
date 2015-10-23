# IMP configuration
We use a config file to pass variables to IMP engine.
The default parameters are visible in `src/config.imp.json`.
You could override some parameters via the file `conf/userconfig.imp.json`.

> Please do not override parameters directly on `src/config.imp.json` as it may be overridden with the next IMP update.

Eventually you could pass a different location for the config file via an environment variable
if you are using snakemake, or via the IMP wrapper script option.



## TODO: document follwing parameters:
{
    "threads": 4,
    "memory_total_gb": 8,
    "memory_per_core_gb": 2,
    "tmp_dir": "tmp",
    "imp_src": "src",
    "raws": {
      "Metagenomics": "",
      "Metatranscriptomics": ""
    },
    "sample": "test",
    "outputdir": "/output",
    "db_path": "/databases",
    "preprocessing_filtering": true,
    "assembler": "idba",
    "trimmomatic": {
        "pkg_url": "https://webdav-r3lab.uni.lu/public/R3lab/IMP/Trimmomatic-Src-0.32.zip",
        "adapter": "TruSeq3",
        "leading": 20,
        "minlen": 40,
        "palindrome_clip_threshold": 30,
        "simple_clip_threshold": 10,
        "trailing": 20,
        "seed_mismatch": 2,
        "window_size": 1,
        "window_quality": 3,
        "strictness": 0.5,
        "target_length": 40,
        "jarfile": "/home/imp/lib/trimmomatic-0.32.jar"
    },
    "idba_ud": {
        "mink": 25,
        "maxk": 99,
        "step": 4,
        "perid": 0.98
    },
    "vizbin": {
        "dimension": 50,
        "kmer": 5,
        "size": 4,
        "theta": 0.5,
        "perp": 30,
        "cutoff": 1000,
        "jarfile": "/home/imp/lib/VizBin-dist.jar"
    },
    "human_filtering": {
        "filter": "hg38",
        "url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    },
    "sortmerna": {
        "pkg_url": "https://webdav-r3lab.uni.lu/public/R3lab/IMP/sortmerna.2.0.tgz",
        "scripts_path": "/home/imp/lib",

        "files": [
            "rfam-5.8s-database-id98",
            "silva-arc-16s-id95",
            "silva-bac-16s-id90",
            "silva-euk-18s-id95",
            "rfam-5s-database-id98",
            "silva-arc-23s-id98",
            "silva-bac-23s-id98",
            "silva-euk-28s-id98"
        ]
    },
    "prokka": {
        "pkg_url": "https://webdav-r3lab.uni.lu/public/R3lab/IMP/prokka-1.11.tar.gz",
        "databases": [
            "cm/Bacteria.i1i",
            "genus/Staphylococcus.phr",
            "hmm/CLUSTERS.hmm.h3f",
            "kingdom/Archaea/sprot.phr"
        ]
    },
    "kegg": {
        "db_ec2pthy": "http://rest.kegg.jp/link/ec/pathway",
        "db_hierarchy": "http://rest.kegg.jp/list/pathway"
    },

}
