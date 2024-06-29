#!/bin/bash

# Make Discinaceae proteomes into a database for BLASTing MAT idiomorphs

# Load modules
module load Bioinformatics
module load blast-plus/2.12.0-fuhtx75

# Set working directory
DB_DIR="$PROJECT_DIR/blastdb"
cd $DB_DIR

# Make BLAST db
for FASTA in *.fasta
do
    echo $FASTA
    # append the file name to the fasta headers so that the BLAST results are informative
    sed "s/>/>${FASTA%.fasta}_/I" $FASTA > ${FASTA%.fasta}_edited.fa
    makeblastdb -in ${FASTA%.fasta}_edited.fa -dbtype prot -out ${FASTA%.fasta}_db
done
blastdb_aliastool -dblist "$(ls *.pto | xargs basename -a -s .pto | xargs)" -dbtype prot -out proteomes_db -title "Discinaceae proteomes DB"

# In case you need to remove all db files generated
# rm *_db.*
