import os
import shutil
import gzip
import json
from copy import deepcopy


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
SRCDIR = os.environ.get("SRCDIR", config['General']['imp_src'])

# get parameters from the command line
DATADIR = os.environ.get("DATADIR", config['General']['datadir'])
OUTPUTDIR = os.environ.get("OUTPUTDIR", config['General']['outputdir'])
MG = os.environ.get("MG", config['General']['raws']['Metagenomics']).split()
MT = os.environ.get("MT", config['General']['raws']['Metatranscriptomics']).split()
SAMPLE = os.environ.get("SAMPLE", config['General']['sample'])
DBPATH = os.environ.get("DBPATH", config['General']['db_path'])
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)

# Get general parameters
THREADS = os.environ.get("THREADS", config['General']['threads'])
MEMTOTAL = os.environ.get("MEMTOTAL", config['General']['memory_total_gb'])
MEMCORE = os.environ.get("MEMCORE", config['General']['memory_per_core_gb'])

# temporary directory will be stored inside the OUTPUTDIR directory
TMPDIR = os.environ.get("TMPDIR", config['General']['tmp_dir'])


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

    if stepname in config and 'pre' in config[stepname]:
        shell(config[stepname]['pre'])
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
