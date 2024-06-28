#!/bin/bash

# Set variables
PROJECT_DIR=$1
SAMPLE_LIST="$PROJECT_DIR/samples.tsv"
OUTPUT_DIR="$PROJECT_DIR/summary"

# Make txt file of categories
printf "Category\n" > $OUTPUT_DIR/categories_antismash.txt

# Run array loop to get list of available categories
while IFS=$'\t' read -r -a files_array
do
    if [ "${files_array[7]}" = "ingroup" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        ANTISMASH_DIR="$PROJECT_DIR/antismash/$NAME"
        cat $ANTISMASH_DIR/index.html | grep '<div class="regbutton ' >> $OUTPUT_DIR/categories_antismash.txt
    fi
done < $LIST

# See available categories
head $OUTPUT_DIR/categories_antismash.txt
cat $OUTPUT_DIR/categories_antismash.txt | sort | uniq
cat $OUTPUT_DIR/categories_antismash.txt | sort | uniq | cut -d " " -f 3,4 | sort | uniq

# Make txt file of counts
printf "Species\tName_full\tName_partial\tNRPS\tPKS\tHybrid\tTerpene\tRiPP\tOther\n" > $OUTPUT_DIR/counts_antismash.txt

# Run array loop
while IFS=$'\t' read -r -a files_array
do
    if [ "${files_array[7]}" = "ingroup" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        ANTISMASH_DIR="$PROJECT_DIR/antismash/$NAME"
        FILE="$ANTISMASH_DIR/index.html"
        SPECIES="${files_array[9]}"
        NAME_FULL="${files_array[15]}"
        NAME_PARTIAL="${files_array[16]}"
        NRPS=$(cat $FILE | grep '<div class="regbutton NRPS ' | wc -l)
        PKS=$(cat $FILE | grep '<div class="regbutton PKS ' | wc -l)
        HYBRID=$(cat $FILE | grep '<div class="regbutton hybrid' | wc -l)
        TERPENE=$(cat $FILE | grep '<div class="regbutton terpene ' | wc -l)
        RIPP=$(cat $FILE | grep '<div class="regbutton RiPP ' | wc -l)
        OTHER=$(cat $FILE | grep '<div class="regbutton other ' | wc -l)
        printf "$LABEL\t$SPECIES\t$NAME_FULL\t$NAME_PARTIAL\t$NRPS\t$PKS\t$HYBRID\t$TERPENE\t$RIPP\t$OTHER\n" >> $OUTPUT_DIR/counts_antismash.txt
    fi
done < $SAMPLE_LIST