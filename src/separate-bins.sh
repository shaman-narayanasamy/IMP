 #!/bin/bash -l

        #!/bin/bash -l
        for cluster in $(tail -n+2 $1 | cut -f2 | sort | uniq)
        do
            # Create directory for the cluster
            mkdir -p Binning/clusterFiles/${cluster}

            # Define headerfiles
            headerfile="Binning/clusterFiles/${cluster}/cluster.${cluster}.contigids"
            protids="Binning/clusterFiles/${cluster}/cluster.${cluster}.protids"

            echo "Obtaining contigs for cluster ${cluster}"
            awk -v c=${cluster} '{{ if($2==c) {{print $1}}  }}' $1 > ${headerfile}
            pullseq -i $2 -n ${headerfile}  > Binning/clusterFiles/${cluster}/cluster.${cluster}.fa

            echo "Obtaining amino acid sequences for cluster ${cluster}"
            awk 'FNR==NR{{a[$1]=$1;next}}{{if(a[$1]) print $4}}' ${headerfile} $3 > ${protids}
            pullseq -i $4 -n ${protids}  > Binning/clusterFiles/${cluster}/cluster.${cluster}.faa 

            echo "Obtaining functional annotations (gff format) for cluster/bin ${cluster}"
            grep -wFf ${headerfile} $5 > Binning/clusterFiles/${cluster}/cluster.${cluster}.gff 
        done
        echo "Complete separating bins"
        touch Binning/separate_bins.done
