#!/bin/bash

# Load modules
module load Bioinformatics
module load sratoolkit/2.10.9-udmejx7

# Get raw reads with SRA Toolkit
mkdir -p $PROJECT_DIR/seqs/raw/sra
cd $PROJECT_DIR/seqs/raw/sra
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ]; then
        
        ACCESSION="${samples_array[5]}"
        prefetch $ACCESSION
        fasterq-dump $ACCESSION
    fi
done < $SAMPLE_LIST

# Move and gzip the files
cd $PROJECT_DIR/seqs/raw
mv $PROJECT_DIR/seqs/raw/sra/*.fastq $PROJECT_DIR/seqs/raw
gzip --fast *.fastq

# Remove the prefetch files
rm -r $PROJECT_DIR/seqs/raw/sra