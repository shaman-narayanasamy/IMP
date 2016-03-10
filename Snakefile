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
if MG and MT:
    include:
        "workflows/integrative/Preprocessing"

    include:
        "workflows/integrative/Assembly"

    include:
        "workflows/integrative/Analysis"

    include:
        "workflows/integrative/Report"

# master command
rule ALL:
    input:
        "preprocessing.done",
        "assembly.done",
        "analysis.done",
        "report.done"
    output:
        touch('workflow.done')
