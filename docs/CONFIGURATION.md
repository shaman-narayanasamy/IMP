# IMP configuration

We use a config file to pass variables to the IMP engine.

The default parameters are visible in `src/config.imp.json`.
You could override some parameters via the user configuration file `conf/userconfig.imp.json`.

> Please do not override parameters directly on `src/config.imp.json` as it may
be overridden with the next IMP update.

Eventually you could pass a different location for the config file via an environment variable
if you are using snakemake, or via the IMP wrapper script option.

~~~
# using snakemake
CONFIGFILE=/path/to/config.json snakemake
# using wrapper script
./IMP -m input/mg.r1.fq -m input/mg.r2.fq -t input/mt.r1.fq -t input/mt.r2.fq -o output_directory -c /path/to/config.json
~~~

Config file must be a valid [JSON](http://json.org) file:

    {
      "threads": 8,
      "output": /home/user/analysis_output_directory,
      "preprocessing_filtering": true
    }

## List of parameters

* threads: Number of max threads to use.
* memory_total_gb: Some tools need to set the max memory they could use.
* memory_per_core_gb: Some tools need to set the max memory they could use per cores.
* tmp_dir: Path to a temporary directory.
* raws - Metagenomics: Path to the metagenomics paired files.
* raws - Metatranscriptomics: Path to the metatranscriptomics paired files.
* outputdir: Path to the output directory.
* db_path: Path to the databases.
* preprocessing_filtering: If you want to filter reads from a database. Can be true or false.
* assembler: The assembler to use. Could be idba or megahit.


### Per tool/step parameters


#### Trimmomatic

* pkg_url: Where to download the trimmomatic package to fetch the adapters databases.
* adapter: What adapter to use.

Following parameters are taken from the [Trimmomatic documentation](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf):

* leading: Cut bases off the start of a read, if below a threshold quality.
* minlen: Specifies the minimum length of reads to be kept.
* palindrome_clip_threshold: Specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
* simple_clip_threshold: Specifies how accurate the match between any adapter etc. sequence must be against a read.
* trailing: Specifies the minimum quality required to keep a base.
* seed_mismatch: specifies the maximum mismatch count which will still allow a full match to be performed.
* window_size: Specifies the number of bases to average across.
* window_quality: Specifies the average quality required.
* strictness: This value, which should be set between 0 and 1, specifies the
balance between preserving as much read length as possible vs. removal of incorrect
bases. A low value of this parameter favours longer reads, while a high value favours read correctness.
* target_length: This specifies the read length which is likely to allow the location of the read within the target sequence to be determined.
* jarfile: Path to the trimmomatic JAR file on your system. (You don't need to set it if you are using the docker container.)


#### idba_ud
* mink: Minimum k value.
* maxk: Maximum k value.
* step: Increment of k-mer of each iteration.
* perid: Similarity for alignment.

#### vizbin
* dimension: 50,
* kmer: 5,
* size: 4,
* theta: 0.5,
* perp: 30,
* cutoff: 1000
* jarfile: Path to the Vizbin JAR file on your system. (You don't need to set it if you are using the docker container.)


#### human_filtering
* filter: Name of the filter.
* url: URL to download database.

#### sortmerna
* pkg_url: Url to download sormerna databases from
* files: Databases to use and index.

#### prokka
* pkg_url: Url to download prokka databases from
* databases: List of databases to use.


#### kegg
* db_ec2pthy and  db_hierarchy: Url to downladod KEgg information from.
