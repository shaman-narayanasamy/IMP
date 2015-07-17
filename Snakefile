include:
    "config"


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


# INCLUDES PROCESSING RULES
include:
    "rules/Util.rules"
include:
    "rules/Preprocessing/master.rules"
include:
    "rules/Assembly/master.rules"
include:
    "rules/Analysis/master.rules"


rule ALL:
    input:
        preprocessing_output_files(),
        assembly_output_files(),
        analysis_output_files()
    shell:
        "echo 'DONE'"
