#!/bin/bash

#SBATCH --job-name=fun_predict
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=150gb
#SBATCH --time=2-00:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

### Memory requirements 
### In my case, using Gyromitra esculenta EST for gene prediction, 50 GB is sufficient for most genomes. 
### For genomes closely related to G. esculenta, more than 150 GB is required. 

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

### Make directory
NAME="$(basename ${GENOME%.fasta})"
mkdir -p $PROJECT_DIR/funannotate/$NAME
cd $PROJECT_DIR/funannotate/$NAME

### Clean and mask assembly
funannotate clean -i $GENOME -o ${NAME}_assembly_cleaned.fasta
funannotate sort -i ${NAME}_assembly_cleaned.fasta -o ${NAME}_assembly_cleaned_sorted.fasta
rm ${NAME}_assembly_cleaned.fasta
funannotate mask -i ${NAME}_assembly_cleaned_sorted.fasta -o ${NAME}_assembly_funmasked.fasta --cpus $CORES

### Predict genes 
###   NOTE: [--busco-db pezizomycotina] argument not recognized by funannotate even though it's in documentation, keep default (don't specify, equals dikarya)
funannotate predict -i ${NAME}_assembly_funmasked.fasta -o fun_out \
--species "$SPECIES" \
--isolate "$ISOLATE" \
--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
--transcript_evidence $PROJECT_DIR/funannotate/Gyresc1_ESTs_20160326_est.fasta \
--cpus $CORES

### Start the submission process with NCBI and get an SBT template file for annotation with antismash, interproscan, and funannotate