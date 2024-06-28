#!/bin/bash

#SBATCH --job-name=busco
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=5gb 
#SBATCH --time=4:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load busco/5.5.0

# Set variables
PROJECT_DIR=$1
SAMPLE=$2
NAME=$(basename ${SAMPLE%.fasta})

# Run BUSCO
cd $PROJECT_DIR/busco
busco -m genome -i "$SAMPLE" -o "$PROJECT_DIR/busco/$NAME" -l $PROJECT_DIR/busco/busco_downloads/lineages/ascomycota_odb10/ --offline