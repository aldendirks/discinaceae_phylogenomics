#!/bin/bash

#SBATCH --job-name=fun_annotate
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=50gb
#SBATCH --time=4:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

### Set environment (run sbatch after activating funannotate conda environment from login terminal)
# conda activate funannotate
source /nfs/turbo/lsa-tyjames/funannotate/envvars.txt # sets paths

### Bash command line input
while getopts p:g:l: flag
do
    case "${flag}" in
        p) PROJECT_DIR=${OPTARG};;
        g) GENOME=${OPTARG};;
        l) LOCUS_TAG=${OPTARG};;
    esac
done
echo "Project directory: $PROJECT_DIR"
echo "FASTA file: $GENOME"
echo "Locus tag prefix: $LOCUS_TAG"

### Set variables and change directory
CORES="$SLURM_NTASKS"
NAME="$(basename ${GENOME%.fasta})"
ANTISMASH_GBK_PATH="$PROJECT_DIR/antismash/$NAME/${NAME}_assembly_funmasked.gbk"
SBT_FILE="$PROJECT_DIR/misc/template.sbt"
cd $PROJECT_DIR/funannotate/$NAME

### Run funannotate (also runs eggnog mapper)
funannotate annotate -i fun_out --antismash $ANTISMASH_GBK_PATH --cpus $CORES --tmpdir $TMPDIR --sbt $SBT_FILE --rename $LOCUS_TAG