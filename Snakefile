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



# Single omics MG workflow
elif MG:
    include:
        "workflows/single_omics/mg/Preprocessing"
    include:
        "workflows/single_omics/mg/Assembly"
    include:
        "workflows/single_omics/mg/Analysis"

    # master command
    rule ALL:
        input:
            "preprocessing.done",
            "assembly.done",
            "analysis.done",
            #"report.done"
        output:
            touch('workflow.done')


elif MT:
    include:
        "workflows/single_omics/mt/Preprocessing"
    include:
        "workflows/single_omics/mt/Assembly"
    include:
        "workflows/single_omics/mt/Analysis"

    # master command
    rule ALL:
        input:
            "preprocessing.done",
            "assembly.done",
            "analysis.done",
            #"report.done"
        output:
            touch('workflow.done')

else:
    raise Exception('No input data.')

# master command
# rule ALL:
#     input:
#         "preprocessing.done",
#         "assembly.done",
#         "analysis.done",
#         "report.done"
#     output:
#         touch('workflow.done')
