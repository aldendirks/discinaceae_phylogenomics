#!/bin/bash

#SBATCH --job-name=mat
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=1gb
#SBATCH --time=1:00:00
#SBATCH --account=alisonhh0
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

### Load modules
echo "$(date +"%T"): Commencing run"
module load Bioinformatics
module load seqkit
module load mafft
module load blast-plus
module load python/3.11.5 # python version I used in the installation of clinker
module load bbtools

### Set variables
CORES="$SLURM_NTASKS"
ANALYSIS_DIR="$PROJECT_DIR/mat"
RESULTS_DIR="$ANALYSIS_DIR/results_$(date +%F)"
BLASTDB="$PROJECT_DIR/blastdb/proteomes_db"
BLAST_QUERY="$PROJECT_DIR/misc/gyromitra-esculenta_mat-locus-16-proteins.faa"
GENBANK_SPLICER_SCRIPT="$PROJECT_DIR/scripts/utilities/slice_genes.py"
GENBANK_EXTRACTOR_SCRIPT="$PROJECT_DIR/scripts/utilities/gbkSPLIT/gbksplit.py"

### Set up directories
mkdir -p $RESULTS_DIR/blast $RESULTS_DIR/phylo $RESULTS_DIR/summary $RESULTS_DIR/clinker/gbk $RESULTS_DIR/coverage
echo "$(date +"%T"): Results directory: $RESULTS_DIR"

### Run BLASTP
echo "$(date +"%T"): Beginning BLASTP"
blastp -db proteomes_db -query $BLAST_QUERY -outfmt 7 -out $RESULTS_DIR/blast/blastp_e-15.out -evalue 1e-15 #-num_threads $CORES
echo "$(date +"%T"): Done BLASTP" # takes about 3 to 4 minutes

### Summarize BLASTP data (presence and absence)

# Split the BLASTP output
cd $RESULTS_DIR/blast
csplit -k $RESULTS_DIR/blast/blastp_e-15.out /BLASTP/ {*} # splits the file into portions between lines containing txt "# BLASTP 2.12.0+"

# Creates xx files; first one, xx00, is empty
rm xx00

# Rename the other xx files
TOTAL=$(ls xx* | tail -n 1 | grep -o "[0-9]*")
for i in $(seq 1 $TOTAL) # use this syntax if you want seq to output leading 0: seq -f %02g 1 $MAX
do
    if [ $i -le 9 ]; then
        FILE="xx0$i"
    else
        FILE="xx$i"
    fi
    mv $FILE blastp_protein$i.out 
done

# Print header of summary file
cd $RESULTS_DIR/summary
printf "Sample\tSuccinate_dehydrogenase\tTranscription_initiation_factor_TFIIE\tFOG_RNA\tCOX13\tAP_endonuclease1\tMAT1-1-1\tHP\tHomeodomain\tUbiquitin\tSHY1\tMAT1-2-1\tUrea_transporter\tHP\tDEAD_box_helicase\tHP\tSLA2\tthallism\n" > $RESULTS_DIR/summary/mat_presence.tsv

# Get presence/absence information for the 16 select genes of the MAT locus
for FASTA in $BLASTDB/*_edited.fa
do
    SAMPLE=$(basename ${FASTA%_edited.fa})
    >$RESULTS_DIR/summary/$SAMPLE.txt
    >$RESULTS_DIR/summary/tmp.txt
    for i in $(seq 1 $TOTAL)
    do
        FILE="$RESULTS_DIR/blast/blastp_protein$i.out"
        cat $FILE | grep $SAMPLE > $RESULTS_DIR/summary/tmp2.txt
        printf "Protein $i\n" >> $RESULTS_DIR/summary/$SAMPLE.txt
        cat $RESULTS_DIR/summary/tmp2.txt >> $RESULTS_DIR/summary/$SAMPLE.txt
        if [ -s tmp2.txt ]; then
            echo "present" >> $RESULTS_DIR/summary/tmp.txt
        else
            echo "absent" >> $RESULTS_DIR/summary/tmp.txt
        fi
    done
    
    # Define thallism
    # Get nth line of file: https://stackoverflow.com/questions/6022384/bash-tool-to-get-nth-line-from-a-file
    MAT1=$(sed '6q;d' $RESULTS_DIR/summary/tmp.txt)
    MAT2=$(sed '11q;d' $RESULTS_DIR/summary/tmp.txt)
    if [ $MAT1 == "absent" ] && [ $MAT2 == "absent" ]; then
        THALLISM="both_absent"
    elif [ $MAT1 == "absent" ] && [ $MAT2 == "present" ]; then
        THALLISM="heterothallic"
    elif [ $MAT1 == "present" ] && [ $MAT2 == "absent" ]; then
        THALLISM="heterothallic"
    elif [ $MAT1 == "present" ] && [ $MAT2 == "present" ]; then
        THALLISM="homothallic"
    else
        echo "something wrong"
    fi
    printf "$SAMPLE\t$(cat $RESULTS_DIR/summary/tmp.txt | tr "\n" "\t")$THALLISM\n" >> $RESULTS_DIR/summary/mat_presence.tsv
done 
rm $RESULTS_DIR/summary/tmp*.txt

# Remove trailing tabs (spaces) in mat_presence.tsv file
# The -i flag tells sed to save the output to the input file itself (in place editing)
# No longer necessary with thallism information
# sed -i 's/[[:space:]]*$//' $RESULTS_DIR/summary/mat_presence.tsv

### Summarize BLASTP data (gene location)

# Print header of summary file
cd $RESULTS_DIR/summary
printf "Sample\tSuccinate_dehydrogenase\tTranscription_initiation_factor_TFIIE\tFOG_RNA\tCOX13\tAP_endonuclease1\tMAT1-1-1\tHP\tHomeodomain\tUbiquitin\tSHY1\tMAT1-2-1\tUrea_transporter\tHP\tDEAD_box_helicase\tHP\tSLA2\tthallism\n" > $RESULTS_DIR/summary/mat_location.tsv

# Get gene location information for the 16 select genes of the MAT locus
for FASTA in $BLASTDB/*_edited.fa
do
    SAMPLE=$(basename ${FASTA%_edited.fa})
    >$RESULTS_DIR/summary/$SAMPLE.txt
    >$RESULTS_DIR/summary/tmp.txt
    for i in $(seq 1 $TOTAL)
    do
        FILE="$RESULTS_DIR/blast/blastp_protein$i.out"
        cat $FILE | grep $SAMPLE > $RESULTS_DIR/summary/tmp2.txt
        printf "Protein $i\n" >> $RESULTS_DIR/summary/$SAMPLE.txt
        cat $RESULTS_DIR/summary/tmp2.txt >> $RESULTS_DIR/summary/$SAMPLE.txt
        if [ $(cat tmp2.txt | wc -l) -gt 0 ]
        then
            cat $RESULTS_DIR/summary/tmp2.txt | cut -f 2 | grep -o "[0-9]*-T1" | sed 's/-T1//' | paste -s -d '_' >> $RESULTS_DIR/summary/tmp.txt
            # GENES=$(cat $RESULTS_DIR/summary/tmp2.txt | cut -f 2 | grep -o "[0-9]*-T1" | sed 's/-T1//' | tr "\n" "," | sed 's/,$//')
            # echo -e "$GENES" >> $RESULTS_DIR/summary/tmp.txt
        else
            echo "NA" >> $RESULTS_DIR/summary/tmp.txt
        fi
    done

    # Define thallism
    # Get nth line of file: https://stackoverflow.com/questions/6022384/bash-tool-to-get-nth-line-from-a-file
    MAT1=$(sed '6q;d' $RESULTS_DIR/summary/tmp.txt)
    echo $MAT1
    MAT2=$(sed '11q;d' $RESULTS_DIR/summary/tmp.txt)
    echo $MAT2
    echo ""
    if [ $MAT1 == "NA" ] && [ $MAT2 == "NA" ]; then
        THALLISM="both_absent"
    elif [ $MAT1 == "NA" ] && [ $MAT2 != "NA" ]; then
        THALLISM="heterothallic"
    elif [ $MAT1 != "NA" ] && [ $MAT2 == "NA" ]; then
        THALLISM="heterothallic"
    elif [ $MAT1 != "NA" ] && [ $MAT2 != "NA" ]; then
        THALLISM="homothallic"
    else
        echo "something wrong"
    fi
    printf "$SAMPLE\t$(cat $RESULTS_DIR/summary/tmp.txt | tr "\n" "\t")$THALLISM\n" >> $RESULTS_DIR/summary/mat_location.tsv
done 
rm $RESULTS_DIR/summary/tmp*.txt

### Identify contiguous gene ranges around MAT alleles - second approach
# This approach gets all ranges instead of just focusing on MAT1-1 and MAT1-2 flanking genes
>$RESULTS_DIR/summary/mat_locus-range_all.tsv
{ 
    read
    while IFS=$'\t' read -ra arrayvar
    do
        >$RESULTS_DIR/summary/tmp9_ranges_gene_list.txt
        IFS=$'\n' arrayvar=($(sort -n <<< "${arrayvar[*]}")); unset IFS # this sorts the array, somehow preserves the sample first
        echo ${arrayvar[0]}
        echo ${arrayvar[@]:1} | tr " " "\n" > $RESULTS_DIR/summary/tmp1.txt
        for GENE in ${arrayvar[@]:1} # this skips the first column (start at column 2, 0 indexing)
        do
            # Define range for that gene
            MIN=$(($GENE - 30))
            MAX=$(($GENE + 30))
            # Get genes that fall within the range
            awk ''$MIN' <= $1 && $1 <= '$MAX'' $RESULTS_DIR/summary/tmp1.txt > $RESULTS_DIR/summary/tmp2.txt
            # Get values of range ends
            MIN_VALUE=$(head -n 1 $RESULTS_DIR/summary/tmp2.txt)
            MAX_VALUE=$(tail -n 1 $RESULTS_DIR/summary/tmp2.txt)
            printf "FUN_$(printf "%06d" $MIN_VALUE):FUN_$(printf "%06d" $MAX_VALUE)\n" >> $RESULTS_DIR/summary/tmp9_ranges_gene_list.txt
        done
        # Get unique values from ranges list
        cat $RESULTS_DIR/summary/tmp9_ranges_gene_list.txt | sort | uniq | tr "\n" "\t" | sed 's/\t$//' > $RESULTS_DIR/summary/tmp10_ranges_gene_list_sorted.txt
        printf "${arrayvar[0]}\t$(cat $RESULTS_DIR/summary/tmp10_ranges_gene_list_sorted.txt)\n" >> $RESULTS_DIR/summary/mat_locus-range_all.tsv
    done 
} < <(cat $RESULTS_DIR/summary/mat_location.tsv | cut -f 1-17 | sed 's/\([0-9][0-9]*\)_\(0\)/\1\t\2/g;s/\t0*/\t/g;s/\tNA/\t/g') # replace commas with tabs, then remove leading 0s, then remove NAs
cat $RESULTS_DIR/summary/mat_locus-range_all.tsv
rm $RESULTS_DIR/summary/tmp*

### Parse Funannotate GenBank file according to range
# Pull genes ranges from GenBank files and output new GenBank files
{
    read
    while IFS=$'\t' read -ra array
    do

        echo -e "\n$(date +"%T"): Beginning extraction of $(echo "${array[0]}")"

        BATCH="$(echo "${array[0]}" | grep -o "Sample_[0-9]*" | grep -o "[0-9]*")"
        SAMPLE="$(echo "${array[0]}" | grep -o "AD-[0-9]*" | grep -o "[0-9]*")"
        echo $BATCH
        echo $SAMPLE

        if [[ $(echo "${array[0]}") == "Sample_8120-AD-10"* ]]
        then
            GBK_FILE="/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/$BATCH-AD_Gyromitra/fastqs_$BATCH-AD/funannotate/Sample_$BATCH-AD-$SAMPLE-filtered/fun_out/predict_results/"*".gbk"
        elif [[ $(echo "${array[0]}") == "Sample_8120-AD-4"* ]]
        then
            # I think the BLAST reference I'm using is from the AJB assembly not the SCGid assembly so later steps are getting messed up... does that make sense?  
            GBK_FILE="/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/$BATCH-AD_Gyromitra/fastqs_$BATCH-AD/funannotate/Sample_$BATCH-AD-$SAMPLE-filtered/fun_out/predict_results/"*".gbk"
        else
            GBK_FILE="/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/$BATCH-AD_Gyromitra/fastqs_$BATCH-AD/funannotate/Sample_$BATCH-AD-$SAMPLE-filtered/fun_out/predict_results/"*".gbk"
        fi

        echo $GBK_FILE

        >$RESULTS_DIR/clinker/gbk_all/$(basename $GBK_FILE)
        for RANGE in "${array[@]:1}"
        do
            if [ $RANGE != "NA" ]; then
                RANGE_WHITESPACE_TRIMMED=$(echo $RANGE | sed 's/ //g')
                echo $RANGE_WHITESPACE_TRIMMED
                python $GENBANK_SPLICER_SCRIPT -i $GBK_FILE -r $RANGE_WHITESPACE_TRIMMED >> $RESULTS_DIR/clinker/gbk_all/$(basename $GBK_FILE)
            fi
        done
    done 
} < $RESULTS_DIR/summary/mat_locus-range_all.tsv

### Run clinker
dos2unix $RESULTS_DIR/clinker/select.txt
mkdir $RESULTS_DIR/clinker/gbk_all_select
while read line
do
    cp $RESULTS_DIR/clinker/gbk_all/$line $RESULTS_DIR/clinker/gbk_all_select
done < $RESULTS_DIR/clinker/select.txt
clinker $RESULTS_DIR/clinker/gbk_all_select/*.gbk -p $RESULTS_DIR/clinker/plot_all_select.html 

### Calculate MAT coverage

# Print header of summary file
cd $RESULTS_DIR/coverage
printf "Full_name\tSpecies\tMAT1-1_gene\tMAT1-1_coverage\tMAT1-2_gene\tMAT1-2-1_coverage\n" > $RESULTS_DIR/coverage/coverage.tsv

# Coverage
{
    read
    while IFS=$'\t' read -ra array
    do
        FULL_NAME="$(echo "${array[0]}")"
        BATCH="$(echo "${array[0]}" | grep -o "Sample_[0-9]*" | grep -o "[0-9]*")"
        SAMPLE="$(echo "${array[0]}" | grep -o "AD-[0-9]*" | grep -o "[0-9]*")"
        BATCH_DIR="/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/$BATCH-AD_Gyromitra/fastqs_$BATCH-AD"

        if [[ $(echo "${array[0]}") == "Sample_8120-AD-10"* ]]
        then
            GBK_FILE="$BATCH_DIR/funannotate/Sample_$BATCH-AD-$SAMPLE-filtered/fun_out/predict_results/"*".gbk"
        elif [[ $(echo "${array[0]}") == "Sample_8120-AD-4"* ]]
        then
            GBK_FILE="$BATCH_DIR/funannotate/Sample_$BATCH-AD-$SAMPLE-filtered/fun_out/predict_results/"*".gbk"
        else
            GBK_FILE="$BATCH_DIR/funannotate/Sample_$BATCH-AD-$SAMPLE-filtered/fun_out/predict_results/"*".gbk"
        fi

        GBK_FILE_REDUCED="$RESULTS_DIR/clinker/gbk/$(basename $GBK_FILE)"
        GBK=$(basename $GBK_FILE)
        SPECIES=${GBK%.gbk}

        echo "$(date +"%T"): Getting FASTA and calculating coverage for $SPECIES"

        # If MAT1-1 exists
        if [[  "${array[6]}" != "NA" ]]
        then

            echo "$(date +"%T"):    Working on MAT1"

            if [[ "${array[6]}" = *"_"* ]]
            then
                echo "${array[6]}" | tr "_" "\n" | sed 's/^/FUN_/' > $RESULTS_DIR/coverage/tmp.txt
            else
                printf "FUN_${array[6]}\n" > $RESULTS_DIR/coverage/tmp.txt
            fi

            $GENBANK_EXTRACTOR_SCRIPT -i $RESULTS_DIR/coverage/tmp.txt -g $GBK_FILE_REDUCED -o $SPECIES-MAT1
            cd $SPECIES-MAT1
            for FASTA1 in *.fasta
            do
                bbmap.sh in=$BATCH_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_1P.fastq.gz \
                in2=$BATCH_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_2P.fastq.gz \
                ref=$FASTA1 \
                nodisk covstats=${FASTA1%.fasta}_covstats.txt \
                covhist=${FASTA1%.fasta}_histogram.txt \
                basecov=${FASTA1%.fasta}_basecov.txt \
                bincov=${FASTA1%.fasta}_bincov.txt 2> ${FASTA1%.fasta}_log.txt
            done

            if [[ "${array[6]}" = *"_"* ]]
            then
                >tmp5.txt
                >tmp6.txt
                for FASTA1 in *.fasta
                do
                    echo ${FASTA1%.fasta} >> tmp5.txt
                    cat ${FASTA1%.fasta}_log.txt | grep "Average coverage:" | grep -o "[0-9]*\.[0-9]*" >> tmp6.txt
                done
                MAT1="$(cat tmp5.txt | tr "\n" "_")"
                MAT1_cov="$(cat tmp6.txt | tr "\n" "_")"
            else
                MAT1="$(echo *.fasta | sed 's/.fasta//')"
                MAT1_cov="$(cat *log.txt | grep "Average coverage:" | grep -o "[0-9]*\.[0-9]*")"
            fi

            cd ..

        else 
            echo "$(date +"%T"):    MAT1 nonexistent"
            MAT1="NA"
            MAT1_cov="NA"
        fi

        # If MAT1-2 exists
        if [[  "${array[11]}" != "NA" ]]
        then

            echo "$(date +"%T"):    Working on MAT2"

            if [[ "${array[11]}" = *"_"* ]]
            then
                echo "${array[11]}" | tr "_" "\n" | sed 's/^/FUN_/' > $RESULTS_DIR/coverage/tmp.txt
            else
                printf "FUN_${array[11]}\n" > $RESULTS_DIR/coverage/tmp.txt
            fi

            $GENBANK_EXTRACTOR_SCRIPT -i $RESULTS_DIR/coverage/tmp.txt -g $GBK_FILE_REDUCED -o $SPECIES-MAT2
            cd $SPECIES-MAT2
            for FASTA2 in *.fasta
            do
                bbmap.sh in=$BATCH_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_1P.fastq.gz \
                in2=$BATCH_DIR/seqs/trimmed/Sample_$BATCH-AD-$SAMPLE-filtered_2P.fastq.gz \
                ref=$FASTA2 \
                nodisk covstats=${FASTA2%.fasta}_covstats.txt \
                covhist=${FASTA2%.fasta}_histogram.txt \
                basecov=${FASTA2%.fasta}_basecov.txt \
                bincov=${FASTA2%.fasta}_bincov.txt 2> ${FASTA2%.fasta}_log.txt
            done

            if [[ "${array[11]}" = *"_"* ]]
            then
                >tmp5.txt
                >tmp6.txt
                for FASTA2 in *.fasta
                do
                    echo ${FASTA2%.fasta} >> tmp5.txt
                    cat ${FASTA2%.fasta}_log.txt | grep "Average coverage:" | grep -o "[0-9]*\.[0-9]*" >> tmp6.txt
                done
                MAT2="$(cat tmp5.txt | tr "\n" "_")"
                MAT2_cov="$(cat tmp6.txt | tr "\n" "_")"
            else
                MAT2="$(echo *.fasta | sed 's/.fasta//')"
                MAT2_cov="$(cat *log.txt | grep "Average coverage:" | grep -o "[0-9]*\.[0-9]*")"
            fi

            cd ..

        else
            echo "$(date +"%T"):    MAT2 nonexistent"
            MAT2="NA"
            MAT2_cov="NA"
        fi

        printf "$FULL_NAME\t$SPECIES\t$MAT1\t$MAT1_cov\t$MAT2\t$MAT2_cov\n" >> $RESULTS_DIR/coverage/coverage.tsv

    done 
} < $RESULTS_DIR/summary/mat_location.tsv
rm $RESULTS_DIR/coverage/tmp*.txt
