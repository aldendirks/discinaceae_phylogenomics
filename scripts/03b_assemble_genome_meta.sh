#!/bin/bash

#SBATCH --job-name=spades_meta
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50gb
#SBATCH --time=6:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load python/3.10.4
module load spades/3.15.5-jhe6qq2

# Set variables
BATCH=$1
SAMPLE=$2

# Run SPAdes meta
#   NOTE: Meta takes error-corrected reads, which was already done. That's why '--only-assembler' is specified. 
#   NOTE: --threads = 16 (default)
mkdir -p "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta"
spades.py --meta --only-assembler -o "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta" \
-1 "$PROJECT_DIR/seqs/corrected_spades/Sample_${BATCH}-AD-${SAMPLE}/corrected/Sample_${BATCH}-AD-${SAMPLE}-filtered_1P.fastq.00.0_0.cor.fastq.gz" \
-2 "$PROJECT_DIR/seqs/corrected_spades/Sample_${BATCH}-AD-${SAMPLE}/corrected/Sample_${BATCH}-AD-${SAMPLE}-filtered_2P.fastq.00.0_0.cor.fastq.gz" \
-s "$PROJECT_DIR/seqs/corrected_spades/Sample_${BATCH}-AD-${SAMPLE}/corrected/Sample_${BATCH}-AD-${SAMPLE}-filtered_1U.fastq.00.0_1.cor.fastq.gz" \
-s "$PROJECT_DIR/seqs/corrected_spades/Sample_${BATCH}-AD-${SAMPLE}/corrected/Sample_${BATCH}-AD-${SAMPLE}-filtered_2U.fastq.00.0_2.cor.fastq.gz"

# Post assembly processing
# Remove contigs smaller than 1000 bp
perl $PROJECT_DIR/scripts/utilities/remove_small_contigs.pl 1000 "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds.fasta" > "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta"
# Run QUAST 
mkdir -p "$PROJECT_DIR/quast/Sample_${BATCH}-AD-${SAMPLE}-meta"
python quast.py -o "$PROJECT_DIR/quast/Sample_${BATCH}-AD-${SAMPLE}-meta" --fungus --rna-finding --split-scaffolds "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta"