#!/bin/bash

#SBATCH --job-name=fastqc_trimmomatic
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=3-00:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Estimated runtime: 1 hour 15 minutes per sample

# Load modules
module load Bioinformatics
module load fastqc/0.11.9-p6ckgle
module load trimmomatic/0.36

# Make directories
mkdir -p $PROJECT_DIR/fastqc/raw $PROJECT_DIR/fastqc/trimmed $PROJECT_DIR/seqs/trimmed

# Run fastqc on raw reads
fastqc -o $PROJECT_DIR/fastqc/raw $PROJECT_DIR/seqs/raw/*

# Run trimmomatic
cd $PROJECT_DIR/seqs/raw
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    ACCESSION="${samples_array[5]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ]; then
        TrimmomaticPE \
        -trimlog $PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE.log \
        -basein $PROJECT_DIR/seqs/raw/${ACCESSION}_1.fastq.gz \
        -baseout $PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered.fastq.gz \
        ILLUMINACLIP:$PROJECT_DIR/misc/TruSeq3-PE-adapters.fa:3:30:10 SLIDINGWINDOW:4:26 LEADING:33 TRAILING:24 MINLEN:50
    fi
done < $SAMPLE_LIST

# Run fastqc on trimmed reads
fastqc -o $PROJECT_DIR/fastqc/trimmed $PROJECT_DIR/seqs/trimmed/*