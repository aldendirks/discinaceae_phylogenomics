#!/bin/bash

# This script recursively makes soft links to proteome files based on data from a tab-delimited text file

# Set destination directory
DST_DIR="$PROJECT_DIR/blastdb"

# Copy proteomes to output dir with sample information
while IFS=$'\t' read -r -a samples_array
do
    if [ "${samples_array[7]}" = "ingroup" ] && [ "${samples_array[0]}" = "1" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        SRC="$PROJECT_DIR/$NAME/fun_out/annotate_results/*proteins.fa"
        echo $SRC
        echo ln -s $SRC $DST_DIR/$NAME.fasta
    fi
done < $SAMPLE_LIST
