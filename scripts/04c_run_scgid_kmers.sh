#!/bin/bash

#SBATCH --job-name=scgid_kmers
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=2gb 
#SBATCH --time=05-00:00:00
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
PROJECT_DIR=$1
BATCH=$2
SAMPLE=$3

# Run scgid kmers module
OUTPUT_DIR="$PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta"
cd "$OUTPUT_DIR"
scgid kmers train -n "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta" --prefix 'meta' --Xmx 2g --mode somoclu --cpus $CORES
scgid kmers annotate -n "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta" --prefix 'meta' --cpus $CORES -s Archaea/Bacteria/Eukaryota^Fungi/Fungi/Viruses