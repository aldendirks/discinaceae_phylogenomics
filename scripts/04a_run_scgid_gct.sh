#!/bin/bash

#SBATCH --job-name=scgid_gct
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=3-00:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load blast-plus/2.12.0-fuhtx75
module load python/3.10.4
module load R/4.2.0

# Set variables
PROJECT_DIR=$1
BATCH=$2
SAMPLE=$3

# Run scgid gct module
OUTPUT_DIR="$PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
scgid gct -n "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta" -g 'Fungi' -sp aspergillus_oryzae --prefix 'meta'