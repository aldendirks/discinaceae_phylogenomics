#!/bin/bash

# Set variables
OUTPUT_DIR="$PROJECT_DIR/summary"
mkdir -p $OUTPUT_DIR/cazy_sample_counts

# Run array loop to get CAZyme counts per sample
while IFS=$'\t' read -r -a files_array
do
    if [ "${files_array[7]}" = "ingroup" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        FILE="$PROJECT_DIR/funannotate/$NAME/fun_out/annotate_misc/annotations.dbCAN.txt"
        cat $FILE | cut -f 3 | sort | uniq -c | sed "s/^[ \t]*//" > $OUTPUT_DIR/cazy_sample_counts/$NAME.txt
    fi
done < $SAMPLE_LIST

# Make txt file
cat $OUTPUT_DIR/cazy_sample_counts/* | cut -d " " -f 2 | sort | uniq > $OUTPUT_DIR/cazy_list.txt
printf "Species\tName_full\tName_partial\t$(cat $OUTPUT_DIR/cazy_list.txt | tr "\n" "\t" | sed 's/[ \t]*$//')\n" > $OUTPUT_DIR/counts_cazy.txt

# Run array loop
while IFS=$'\t' read -r -a files_array
do
    if [ "${files_array[7]}" = "ingroup" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        SPECIES="${files_array[9]}"
        NAME_FULL="${files_array[15]}"
        NAME_PARTIAL="${files_array[16]}"
        >$OUTPUT_DIR/tmp.txt
        while read line
        do
            COUNT=$(cat $OUTPUT_DIR/cazy_sample_counts/$NAME.txt | grep -w $line | grep -o "[0-9]* ")
            if [ -z "$COUNT" ]; then
                echo "0" >> $OUTPUT_DIR/tmp.txt
            else
                echo $COUNT >> $OUTPUT_DIR/tmp.txt
            fi
        done < $OUTPUT_DIR/cazy_list.txt
        printf "$SPECIES\t$NAME_FULL\t$NAME_PARTIAL\t$(cat $OUTPUT_DIR/tmp.txt | tr "\n" "\t" | sed 's/[ \t]*$//')\n" >> $OUTPUT_DIR/counts_cazy.txt
    fi
done < $SAMPLE_LIST

rm $OUTPUT_DIR/tmp.txt