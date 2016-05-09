#! /bin/bash -l

# Set the environment
#module use /home/users/claczny/.local/resif/mysoftware/modules/all/ #TODO check if needed
#module use /work/projects/ecosystem_biology/local_tools/privatemodules/
#module use /opt/apps/resif/devel/v1.0-20150402/lcsb/all/
#module use $RESIF_ROOTINSTALL/lcsb/modules/all/

#module load AMPHORA2/parallelized3
#module load AMPHORA2
# Check for command line args.

AMPHORA2_home=/home/imp/lib/AMPHORA2/

OPTIND=1 #reset for getopts
INFILE=""
CPUS=1
OUTDIR=""
DNA="-DNA"

function show_help { echo "Usage: `basename $0` -i <INFILE (fasta)> -o <OUTSUMMARY (tsv)> [mandatory] -d <OUTDIR> -c <CPUS> -p [optional]" >&2 ; exit 1; }

while getopts "h:?:d:c:i:o:p" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    d)  OUTDIR=${OPTARG}
        ;;
    c)  CPUS=${OPTARG}
        ;;
    i)  INFILE=${OPTARG}
        ;;
    o)  OUTSV=${OPTARG}
        ;;
    p)  DNA=""
	;;
    esac
done

if [ -z "${INFILE}" ] || [ -z "${CPUS}" ] || [ -z "${OUTSV}" ];
then
  show_help
fi

CPUSTRING=""
if [ $CPUS -gt 1 ];
then
  CPUSTRING="-CPUs ${CPUS}"
fi

#run commands: MarkerScanner MarkAlignTrim Phylotyping for fasta with bacterial dna, contig or genome
OLDIR=$PWD

if [ -z "$OUTDIR" ];
then
  CMD="perl $AMPHORA2_home/Scripts/MarkerScanner.pl ${DNA} -Bacteria ${INFILE}"
  echo ${CMD}
  eval ${CMD}
  CMD="perl $AMPHORA2_home/Scripts/MarkerAlignTrim.pl -WithReference ${CPUSTRING}"
  echo ${CMD}
  eval ${CMD}
  CMD="perl $AMPHORA2_home/Scripts/Phylotyping.pl ${CPUSTRING} > ${OLDIR}/${OUTSV}"
  echo ${CMD}
  eval ${CMD}
else
  CMD="perl $AMPHORA2_home/Scripts/MarkerScanner.pl ${DNA} -Bacteria -Outdir ${OUTDIR} ${INFILE}"
  echo ${CMD}
  eval ${CMD}
  CMD="perl $AMPHORA2_home/Scripts/MarkerAlignTrim.pl -WithReference ${CPUSTRING} -Directory ${OUTDIR}"
  echo ${CMD}
  eval ${CMD}
  echo "changing dir to ${OUTDIR}"
  cd ${OUTDIR}
  CMD="perl $AMPHORA2_home/Scripts/Phylotyping.pl ${CPUSTRING} > ${OUTSV}"
  echo ${CMD}
  eval ${CMD}
  echo "changing dir to ${OLDIR}"
  cd ${OLDIR}
fi

