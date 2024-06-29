#!/bin/bash

# Set variables
OUTPUT_DIR="$PROJECT_DIR/summary"

# Make txt file
printf "Species\tName_full\tName_partial\tSize\tGC\tGenes\tAverage_gene_length\ttRNA\tPFAM\tMEROPS\tSecretion\tAverage_exon_length\tAvg_protein_length\tExons\tTransposable_elements\n" > $OUTPUT_DIR/genome_stats.txt

# Run array loop to get various genome stats
while IFS=$'\t' read -r -a files_array
do
    if [ "${files_array[7]}" = "ingroup" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        FUN_DIR="$PROJECT_DIR/funannotate/$NAME/fun_out/annotate_results"
        FUN_DIR_TE="$PROJECT_DIR/funannotate/$NAME/fun_out/logfiles/funannotate-predict.log"
        SPECIES="${files_array[9]}"
        NAME_FULL="${files_array[15]}"
        NAME_PARTIAL="${files_array[16]}"
        SIZE="$(cat $FUN_DIR/*.stats.json | grep '"length":' | grep -o "[0-9]*")"
        GC="$(cat $FUN_DIR/*.stats.json | grep '"GC_content":' | grep -o "[0-9]*\\.[0-9]*")"
        GENES="$(cat $FUN_DIR/*.stats.json | grep '"genes":' | grep -o "[0-9]*")"
        AVG_GENE_LENGTH="$(cat $FUN_DIR/*.stats.json | grep '"avg_gene_length":' | grep -o "[0-9]*\\.[0-9]*")"
        TRNA="$(cat $FUN_DIR/*.stats.json | grep '"tRNA":' | grep -o "[0-9]*")"
        PFAM="$(cat $FUN_DIR/*.stats.json | grep '"pfam":' | grep -o "[0-9]*")"
        MEROPS="$(cat $FUN_DIR/*.stats.json | grep '"merops":' | grep -o "[0-9]*")"
        SECRETION="$(cat $FUN_DIR/*.stats.json | grep '"secretion":' | grep -o "[0-9]*")"
        AVG_EXON_LENGTH="$(cat $FUN_DIR/*.stats.json | grep '"avg_exon_length":' | grep -o "[0-9]*\\.[0-9]*")"
        AVG_PROTEIN_LENGTH="$(cat $FUN_DIR/*.stats.json | grep '"avg_protein_length":' | grep -o "[0-9]*\\.[0-9]*")"
        EXONS="$(cat $FUN_DIR/*.stats.json | grep '"total_exons":' | grep -o "[0-9]*")"
        TE="$(cat $FUN_DIR_TE | grep "transposable elements" | tr -d "," | grep -o "; [0-9]* transposable elements" | grep -o "[0-9]*")"
        printf "$SPECIES\t$NAME_FULL\t$NAME_PARTIAL\t$SIZE\t$GC\t$GENES\t$AVG_GENE_LENGTH\t$TRNA\t$PFAM\t$MEROPS\t$SECRETION\t$AVG_EXON_LENGTH\t$AVG_PROTEIN_LENGTH\t$EXONS\t$TE\n" >> $OUTPUT_DIR/genome_stats.txt
    fi
done < $SAMPLE_LIST