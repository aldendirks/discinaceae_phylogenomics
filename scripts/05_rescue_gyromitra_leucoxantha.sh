#!/bin/bash

#SBATCH --job-name=rescue_gyromitra_leucoxantha
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50gb
#SBATCH --time=12:00:00
#SBATCH --account=alisonhh0
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load blast-plus/2.12.0-fuhtx75
module load python/3.10.4
module load seqtk/1.3
module load spades/3.15.5-jhe6qq2

# Set variables
CORES="$SLURM_NTASKS"
PROJECT_DIR=$1
BATCH="6977"
SAMPLE="22"
REFERENCE="$PROJECT_DIR/scgid/Sample_4981-AD-7-scgid-meta/meta_scgid_output/consensus/meta.consensus.filtered.assembly.fasta"
TRIMMED_FASTQ="$PROJECT_DIR/seqs/trimmed/Sample_6977-AD-22-filtered"*
OUTPUT_DIR="$PROJECT_DIR/rescue"
DB_DIR="$OUTPUT_DIR/blastdb"
BLAST_DIR="$OUTPUT_DIR/blast-out"
FASTQ_DIR="$OUTPUT_DIR/fastq"
mkdir -p $DB_DIR $BLAST_DIR $FASTQ_DIR

# Make BLAST database for Gyromitra leucoxantha MICH352087
ln -s $REFERENCE $DB_DIR/Gyromitra-leucoxantha_ACD0366.fasta
makeblastdb -in $DB_DIR/Gyromitra-leucoxantha_ACD0366.fasta -dbtype nucl -out $DB_DIR/Gyromitra-leucoxantha_ACD0366_db

# BLAST reads and filter
for FASTQ in $TRIMMED_FASTQ
do
    # Get file basename without file extension
    NAME="$(basename ${FASTQ%.fastq.gz})"

    # Convert fastq to fasta
    seqtk seq -a $FASTQ > $BLAST_DIR/tmp.fasta 

    # BLAST the trimmed reads against the reference genome
    blastn -query $BLAST_DIR/tmp.fasta -db $DB_DIR/Gyromitra-leucoxantha_ACD0366_db -outfmt '6' -evalue 1e-5 \
    -out $BLAST_DIR/${NAME}_blast.out -num_threads $CORES

    # Get read names
    cat $BLAST_DIR/${NAME}_blast.out | cut -f 1 > $BLAST_DIR/${NAME}_blast.txt

    # Subset FASTQ
    seqtk subseq $FASTQ $BLAST_DIR/${NAME}_blast.txt > $FASTQ_DIR/${NAME}_subset.fastq

    # Remove temporary file
    rm $BLAST_DIR/tmp.fasta
done

# Match up paired reads
# Make temporary file
> $FASTQ_DIR/tmp.txt
# Add all fasta headers to temporary file
for FASTQ in $FASTQ_DIR/Sample_6977-AD-22-filtered_1P_subset.fastq $FASTQ_DIR/Sample_6977-AD-22-filtered_2P_subset.fastq
do
    cat $FASTQ | grep "@" | cut -d " " -f 1 | sed 's/@//' | sort | uniq >> $FASTQ_DIR/tmp.txt
done
# Get list of unique headers
cat $FASTQ_DIR/tmp.txt | sort | uniq > $FASTQ_DIR/all_reads.txt
# Reextract sequences from FASTQ files 
seqtk subseq $(echo $TRIMMED_FASTQ | tr " " "\n" | grep "1P") $FASTQ_DIR/all_reads.txt > $FASTQ_DIR/Sample_6977-AD-22-filtered_1P_matched.fastq
seqtk subseq $(echo $TRIMMED_FASTQ | tr " " "\n" | grep "2P") $FASTQ_DIR/all_reads.txt > $FASTQ_DIR/Sample_6977-AD-22-filtered_2P_matched.fastq
# Remove files
rm $FASTQ_DIR/tmp.txt $FASTQ_DIR/all_reads.txt

# gzip files
gzip --fast $FASTQ_DIR/*.fastq

# Run SPAdes
mv $PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta $PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-fungus-contaminated
mkdir -p "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta"
spades.py --meta --only-assembler -o "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta" \
-1 "$FASTQ_DIR/Sample_6977-AD-22-filtered_1P_matched.fastq.gz" \
-2 "$FASTQ_DIR/Sample_6977-AD-22-filtered_2P_matched.fastq.gz" \
-s "$FASTQ_DIR/Sample_6977-AD-22-filtered_1U.fastq.gz" \
-s "$FASTQ_DIR/Sample_6977-AD-22-filtered_2U.fastq.gz"

# Post assembly processing
# Remove contigs smaller than 1000 bp
perl $PROJECT_DIR/misc/remove_small_contigs.pl 1000 "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds.fasta" > "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta"
# Run QUAST 
mv $PROJECT_DIR/quast/Sample_${BATCH}-AD-${SAMPLE}-meta $PROJECT_DIR/quast/Sample_${BATCH}-AD-${SAMPLE}-fungus-contaminated
mkdir -p "$PROJECT_DIR/quast/Sample_${BATCH}-AD-${SAMPLE}-meta"
python quast.py -o "$PROJECT_DIR/quast/Sample_${BATCH}-AD-${SAMPLE}-meta" --fungus --rna-finding --split-scaffolds "$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta"

# Finally, modify the SCGid output so the pipeline can find the right genome
mv $PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta $PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-fungus-contaminated
mkdir -p "$PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta/meta_scgid_output/consensus"
cp $PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-meta/scaffolds_trimmed.fasta $PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta/meta_scgid_output/consensus/meta.consensus.filtered.assembly.fasta