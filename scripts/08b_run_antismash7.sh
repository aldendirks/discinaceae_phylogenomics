#!/bin/bash

#SBATCH --job-name=antismash
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=20gb
#SBATCH --time=4:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

### Make sure to load conda environment before running this script
# conda activate antismash

### Show antismash commands, version 7.1.0
# antismash --help-showall

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

### Set variables and change directory
CORES="$SLURM_NTASKS"
NAME="$(basename ${GENOME%.fasta})"
mkdir -p $PROJECT_DIR/antismash/$NAME
cd $PROJECT_DIR/antismash/$NAME
GFF3_FILE="$PROJECT_DIR/funannotate/$NAME/fun_out/predict_results/*.gff3"
FASTA_FILE="$PROJECT_DIR/funannotate/$NAME/*funmasked.fasta"

### Run antismash
antismash --taxon fungi --cassis --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --cb-subclusters --cb-knownclusters --pfam2go --rre --smcog-trees --tfbs \
--output-dir . --cpus $CORES --genefinding-gff3 $GFF3_FILE $FASTA_FILE

### Convert GenBank output files to protein FASTA files
for FILE in *.gbk
do
    HEADER="$SPECIES $ISOLATE ${FILE%.region00*.gbk}"
    HEADER="${HEADER// /_}"
    python $PROJECT_DIR/scripts/utilities/gbk_converter.py $FILE $HEADER
done
FILE_NAME="$SPECIES $ISOLATE"
FILE_NAME="${FILE_NAME// /_}"
cat *region*.fasta > "$FILE_NAME"_antismash.fasta