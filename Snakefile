# include configuration file
include:
    "rules/ini/config"

# define the data types used and the assembly
if MG and MT:
    TYPES = ['mg', 'mt']
    ASS = 'mgmt'
elif MG:
    TYPES = ['mg']
    ASS = 'mg'
elif MT:
    TYPES = ['mt']
    ASS = 'mt'

workdir:
    OUTPUTDIR

# include rules for the workflow based on the input parameters


# INTEGRATIVE MG-MT workflow
if MG and MT:
    include:
        "workflows/integrative/Preprocessing"

    include:
        "workflows/integrative/Assembly"

    include:
        "workflows/integrative/Analysis"

    include:
        "workflows/integrative/Binning"

    include:
        "workflows/integrative/Report"


# Single omics MG workflow
elif MG:
    include:
        "workflows/single_omics/mg/Preprocessing"
    include:
        "workflows/single_omics/mg/Assembly"
    include:
        "workflows/single_omics/mg/Analysis"
    include:
        "workflows/single_omics/mg/Binning"
    include:
        "workflows/single_omics/mg/Report"


elif MT:
    include:
        "workflows/single_omics/mt/Preprocessing"
    include:
        "workflows/single_omics/mt/Assembly"
    include:
        "workflows/single_omics/mt/Analysis"
    include:
        "workflows/single_omics/mt/Binning"
    include:
        "workflows/single_omics/mt/Report"

else:
    raise Exception('No input data.')

# master command
rule ALL:
    input:
        "preprocessing.done",
        "assembly.done",
        "analysis.done",
        "report.done",
        "binning.done"
    output:
        touch('workflow.done')
