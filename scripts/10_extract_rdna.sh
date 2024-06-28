#!/bin/bash

#SBATCH --job-name=extractITS
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=10gb 
#SBATCH --time=10:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load bedtools2/2.30.0-svcfwbm
module load blast-plus/2.12.0-fuhtx75
module load python/3.10.4

# Set variables
PROJECT_DIR=$1
SAMPLE_LIST="$PROJECT_DIR/samples.tsv"

# Extract rDNA
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    GENOME="$PROJECT_DIR/genomes/$NAME.fasta"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ]; then

        # Make output directory and run extractITS.py
        OUTPUT_DIR="$PROJECT_DIR/rdna/Sample_${BATCH}-AD-${SAMPLE}"
        mkdir -p $OUTPUT_DIR
        python extractITS.py -i $GENOME -o $OUTPUT_DIR -which all -name $NAME

        # Rename FASTA files
        for FILE in $OUTPUT_DIR/$NAME.SSU.fasta $OUTPUT_DIR/$NAME.full.fasta $OUTPUT_DIR/$NAME.LSU.fasta
        do 
            sed -i "s/>.*/>$NAME/" > $FILE
        done 
    fi
done < $SAMPLE_LIST