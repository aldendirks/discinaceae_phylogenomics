#!/bin/bash

#SBATCH --job-name=scgid_kmers-extract_consensus
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1gb 
#SBATCH --time=1:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load blast-plus/2.12.0-fuhtx75
module load cuda/11.5.1
module load openmpi/4.1.4
module load python/3.10.4
module load R/4.2.0

# Set variables
CORES="$SLURM_NTASKS"
BATCH=$1
SAMPLE=$2
CLASS_ID=6 # you will need to change this number depending on the class given to your manual ESOM selection. 

# Run scgid kmers extract module and get consensus
OUTPUT_DIR="$PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta"
cd "$OUTPUT_DIR"
scgid kmers extract -n "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta" --prefix 'meta' -c meta_scgid_output/kmers/selection.cls -cid $CLASS_ID
scgid consensus -n "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta" --prefix 'meta' --venn --exclude_annot_nt