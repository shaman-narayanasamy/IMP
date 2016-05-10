#! /bin/bash -l

###
#   Script to prepare proteins for downstream analysis using phylophlan-impute
#   Cedric C. Laczny (2015)
###

EXPECTED_ARGS=2
E_BADARGS=65

if (( $# != $EXPECTED_ARGS ))
then
    echo "Usage: `basename $0` <directory w/ unprepared .faa files> <OUT_DIR>"
    exit $E_BADARGS
fi

IN_DIR="${1}"
OUT_DIR="${2}"

date

echo "STATUS: Preparing phylophlan-impute input ..."

#create outdir unless exists
if [[ -d ${OUT_DIR} ]]; then
    echo "INFO: DEST_DIR (${OUT_DIR}) already exists."
else
    mkdir -p ${OUT_DIR}
    echo "INFO: Created ${OUT_DIR}."
fi

for i in $(find ${IN_DIR} -name "*.faa" -maxdepth 1); do
   CMD="cp -a ${i} ${OUT_DIR}"
   echo "${CMD}"
   eval "time ${CMD}"
done

echo "STATUS: Replacing whitespaces with underscores ..."
CMD="for j in ${OUT_DIR}/*.faa; do sed -i \"/^>/{s/ /_/g}\" \"\$j\"; done"
echo "${CMD}"
eval "time ${CMD}"
echo "done (replacing whitespaces)."

echo "STATUS: Prepending filenames to FASTA IDs ..."
CMD="for j in ${OUT_DIR}/*.faa; do NAME=\$(basename \"\$j\" .faa); sed -i \"s/^>/>\${NAME}_/\" \"\$j\"; done"
echo "${CMD}"
eval "time ${CMD}"
echo "done (prepending filename)."

echo "STATUS: Removing the asterisk(s) in and at the ends of the aminoacid sequences ..."
# Remove the asterisk at the end of a peptide -> PhyloPhlAn complains otherwise. Also remove * within sequence
CMD="for j in ${OUT_DIR}/*.faa; do sed -i 's/\*//g' \"\$j\"; done"
echo "${CMD}"
eval "time ${CMD}"
echo "done (removing asterisk)."

#touch prepareProteinsForPhylophlanImpute.done

echo "done (preparing phylophlan-impute input)."

date

