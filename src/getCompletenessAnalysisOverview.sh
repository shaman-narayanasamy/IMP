#! /bin/bash -l

###
#   Script to extract and summarize the per-CAG essential gene analysis results
#   Cedric C. Laczny (2015)
###
#### Note that this needs per cluster .hits file

EXPECTED_ARGS=2
E_BADARGS=65

if (( $# != $EXPECTED_ARGS ))
then
    echo "Usage: `basename $0` <directory w/ essential.hits files> <OUT_FILE (csv)>"
    exit $E_BADARGS
fi

IN_DIR="${1}"
OUT_FILE="${2}"
COUNTER=0

date

echo "STATUS: Collecting completeness results..."
echo "#ClusterName,num. uniques,num. multiples" > ${OUT_FILE}
for IN_FILE in $(find ${IN_DIR} -name "*essential.hits" | sort)
do
    COUNTER=$(($COUNTER + 1))
    (echo -ne "$(basename "${IN_FILE}" .essential.hits)," && grep -v "^#" ${IN_FILE} | awk '{print $4}' | sort | uniq -c | awk 'BEGIN { uniques=0; multis=0 } {if($1 == 1) uniques += 1; else multis +=1;} END { print uniques "," multis }') >> ${OUT_FILE}
done

echo "done (collecting ${COUNTER} completeness results)."

date
