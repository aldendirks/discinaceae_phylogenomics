#!/bin/bash

#SBATCH --job-name=correct_reads
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=100gb
#SBATCH --time=4:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Estimated runtime: 1 hour 40 minutes per sample

# Load modules
module load Bioinformatics
module load spades/3.15.5-jhe6qq2

# Set variables
PROJECT_DIR=$1
BATCH=$2
SAMPLE=$3

# Run SPAdes error correction
mkdir -p $PROJECT_DIR/seqs/corrected_spades/Sample_$BATCH-AD-$SAMPLE
spades.py --only-error-correction -o $PROJECT_DIR/seqs/corrected_spades/Sample_$BATCH-AD-$SAMPLE \
-1 $PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_1P.fastq.gz \
-2 $PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_2P.fastq.gz \
-s $PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_1U.fastq.gz \
-s $PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_2U.fastq.gz