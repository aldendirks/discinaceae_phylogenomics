#!/bin/bash

#SBATCH --job-name=interproscan
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=100gb
#SBATCH --time=8:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

### Set environment (run sbatch after activating funannotate conda environment from login terminal)
# conda activate funannotate
source /nfs/turbo/lsa-tyjames/funannotate/envvars.txt # sets paths

### Set environment (run sbatch after activating funannotate conda environment from login terminal)
# conda activate funannotate
source /nfs/turbo/lsa-tyjames/funannotate/envvars.txt # sets paths

### Bash command line input
while getopts s:i:g: flag
do
    case "${flag}" in
        s) SPECIES=${OPTARG};;
        i) ISOLATE=${OPTARG};;
        g) GENOME=${OPTARG};;
    esac
done
echo "Species: $SPECIES"
echo "Isolate/accession: $ISOLATE"
echo "FASTA file: $GENOME"

### Set variables
CORES="$SLURM_NTASKS"

### Run interproscan
NAME="$(basename ${GENOME%.fasta})"
mkdir -p $PROJECT_DIR/funannotate/$NAME
cd $PROJECT_DIR/funannotate/$NAME
funannotate iprscan -i fun_out -m local --cpus $(($CORES-2))
