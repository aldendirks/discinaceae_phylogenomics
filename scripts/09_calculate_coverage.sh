#!/bin/bash

#SBATCH --job-name=bbmap
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb 
#SBATCH --time=48:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load bbtools/38.96

# Get coverage
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ]; then
        OUTPUT_DIR="$PROJECT_DIR/coverage/Sample_$BATCH-AD-$SAMPLE"
        mkdir -p $OUTPUT_DIR
        cd $OUTPUT_DIR
        bbmap.sh in=$PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_1P.fastq.gz \
        in2=$PROJECT_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_2P.fastq.gz \
        ref=$PROJECT_DIR/genomes/$NAME.fasta \
        nodisk covstats=covstats.txt covhist=histogram.txt basecov=basecov.txt bincov=bincov.txt 2> log.txt
    fi
done < $SAMPLE_LIST