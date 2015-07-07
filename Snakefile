import os
import shutil
import gzip
import json
import bz2
from copy import deepcopy
import subprocess


def dict_merge(a, b):
    """
    Deep merge 2 dicts together
    """
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.items():
        if k in result and isinstance(result[k], dict):
            result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result

# default configuration file
configfile:
    "src/config.imp.json"

# default executable for snakmake
shell.executable("bash")

# custom configuration file
CUSTOM_CONFIG_PATH = "conf/userconfig.imp.json"

# merge 2 configurations files together
if os.path.exists(CUSTOM_CONFIG_PATH):
    with open(CUSTOM_CONFIG_PATH, 'r') as rhandle:
        data = json.load(rhandle)
        config = dict_merge(config, data)


# some parameters
SRCDIR = os.environ.get("SRCDIR", config['imp_src'])

KOOPA = os.environ.get("KOOPA", None)
# get parameters from the command line
OUTPUTDIR = os.environ.get("OUTPUTDIR", config['outputdir'])
MG = os.environ.get("MG", config['raws']['Metagenomics']).split()
MT = os.environ.get("MT", config['raws']['Metatranscriptomics']).split()
SAMPLE = os.environ.get("SAMPLE", config['sample'])
DBPATH = os.environ.get("DBPATH", config['db_path'])
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)

# Get general parameters
THREADS = os.environ.get("THREADS", config['threads'])
MEMTOTAL = os.environ.get("MEMTOTAL", config['memory_total_gb'])
MEMCORE = os.environ.get("MEMCORE", config['memory_per_core_gb'])

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = os.environ.get("TMPDIR", config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)


def prepare_environment(stepname):
    """
    Prepare the output directories and logs.
    stepname: the name of the pipeline step
    return: the step master directory, the step log
    """
    out = os.path.join(OUTPUTDIR, stepname)
    # mkdirs
    if not os.path.exists(out):
        os.makedirs(out)
    elif not os.path.isdir(out):
        raise OSError("//[IMP] Output is not a directory: %s" % out)
    if not os.path.exists(TMPDIR):
        os.makedirs(TMPDIR)
    bench = os.path.join(out, 'benchmarks')
    if not os.path.exists(bench):
        os.makedirs(bench)

    return out, os.path.join(out, '%s.log' % stepname)


# INCLUDES RULES
include:
    "rules/Util.rules"
include:
    "rules/Preprocessing/master.rules"
include:
    "rules/Assembly/master.rules"
include:
    "rules/Analysis/master.rules"


# locate source directory and name scripts
src = lambda p: os.path.join(SRCDIR, p)


rule ALL:
    input:
        preprocessing_output_files(),
        assembly_output_files(),
        analysis_output_files()

    shell:
        "echo 'DONE'"

rule MODULE_LOAD_TEST:
    shell:
        """
        IMPPRL="{config[preload][test]}"; if [[ -n $IMPPRL ]]; then $IMPPRL; fi
        """

if KOOPA:
    print(KOOPA)

rule T:
    shell:
        """
        echo {MT}
        """
