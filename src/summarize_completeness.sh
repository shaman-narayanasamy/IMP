
        #src/getCompletenessAnalysisOverview.sh Binning/clusterFiles Binning/CompletenessSummary.csv
        IN_DIR="Binning/clusterFiles"
        OUT_FILE=Binning/CompletenessSummary.csv
        COUNTER=0
        
        echo "STATUS: Collecting completeness results..."
        echo "#ClusterName,num. uniques,num. multiples" > ${OUT_FILE}
        for IN_FILE in $(find $IN_DIR -maxdepth 1 -name "*essential.hits" | sort)
        do
            COUNTER=$(($COUNTER + 1))
            (echo -ne "$(basename "$IN_FILE" .essential.hits)," && grep -v "^#" $IN_FILE | awk '{print $4}' | sort | uniq -c | awk 'BEGIN { uniques=0; multis=0 } {if($1 ==1) uniques += 1; else multis +=1;} END { print uniques "," multis }') >> ${OUT_FILE}
        done
        
        echo "done (collecting ${COUNTER} completeness results)."
