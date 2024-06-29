#!/bin/bash

#SBATCH --job-name=phylogenomics
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=100gb 
#SBATCH --time=5-00:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load python/3.10.4
module load muscle/3.8.31
module load fasttree/2.1.10

# Set variables
CORES="$SLURM_NTASKS"
OUTPUT_DIR="$PROJECT_DIR/phylogenomics"
RUNS_DIR="$OUTPUT_DIR/runs"
BUSCO_PHYLO_RESULTS_DIR="$OUTPUT_DIR/busco_phylogenomics_results"
IQTREE_RESULTS_DIR="$OUTPUT_DIR/iqtree_results"
mkdir -p $RUNS_DIR/all $RUNS_DIR/select $BUSCO_PHYLO_RESULTS_DIR $IQTREE_RESULTS_DIR

# Soft link BUSCO runs - all
while IFS=$'\t' read -r -a files_array
do
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    RUN="$PROJECT_DIR/busco/$NAME/run_ascomycota_odb10"
    ln -s $RUN $RUNS_DIR/all/run_$NAME
done < $SAMPLE_LIST

# Soft link BUSCO runs - select
while IFS=$'\t' read -r -a files_array
do
    if [ "${files_array[0]}" = "1" ]; then
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        RUN="$PROJECT_DIR/busco/$NAME/run_ascomycota_odb10"
        ln -s $RUN $RUNS_DIR/select/run_$NAME
    fi
done < $SAMPLE_LIST

# Run BUSCO phylogenomics pipeline
for RUNS in $RUNS_DIR/all $RUNS_DIR/select
do
    # concatenation
    # --stop_early generates the files that are used by IQTREE2 without continuing on to phylogenetic inference
    python BUSCO_phylogenomics.py -d $RUNS -o $BUSCO_PHYLO_RESULTS_DIR/supermatrix_$(basename $RUNS)_90 --supermatrix --threads $CORES --percent_single_copy 90 --stop_early
    python BUSCO_phylogenomics.py -d $RUNS -o $BUSCO_PHYLO_RESULTS_DIR/supermatrix_$(basename $RUNS)_75 --supermatrix --threads $CORES --percent_single_copy 75 --stop_early
    # coalescence
    python BUSCO_phylogenomics.py -d $RUNS -o $BUSCO_PHYLO_RESULTS_DIR/supertree_$(basename $RUNS) --supertree --threads $CORES
    # this exact formatting (no spaces between -D and quotation marks, space after =, etc.) was found to be necessary for some reason
    java -D"java.library.path= "$ASTRAL_LIB_PATH"" -jar $ASTRAL_JAR_PATH \
    -i $BUSCO_PHYLO_RESULTS_DIR/supertree_$(basename $RUNS)/ALL.trees -o $BUSCO_PHYLO_RESULTS_DIR/supertree_$(basename $RUNS)/astral.tree
done

# Run IQTREE2 pipeline
for RUNS in $RUNS_DIR/all $RUNS_DIR/select
do
    for CUTOFF in 75 90
    do
        # Set up directories
        IQTREE_ANALYSIS_DIR="$IQTREE_RESULTS_DIR/$(basename $RUNS)_$CUTOFF"
        mkdir -p $ANALYSIS_DIR
        SUPERMATRIX_ALN_DIR="$BUSCO_PHYLO_RESULTS_DIR/supermatrix_$(basename $RUNS)_$CUTOFF/trimmed_alignments"
        ALN_DIR="$ANALYSIS_DIR/trimmed_alignments"
        ln -s $SUPERMATRIX_ALN_DIR $ALN_DIR
        cd $IQTREE_ANALYSIS_DIR

        # Infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
        iqtree2 -p $ALN_DIR -B 1000 -mset WAG,LG,JTT --prefix concat -T $CORES

        # Test phylogenetic assumptions
        #   stationarity: amino acid frequencies remain constant over time
        #   homogeneity: substitution rates remain constant over time
        iqtree2 -p $ALN_DIR --symtest-remove-bad -B 1000 -mset WAG,LG,JTT --prefix good_partitions -T $CORES

        # Compare trees from all partitions and from "good" partitions only
        cat concat.treefile good_partitions.treefile > species_trees.treefile
        iqtree2 -p $ALN_DIR -z species_trees.treefile -n 0 -zb 10000 -zw -au -mset WAG,LG,JTT --prefix topotest_symtest -T $CORES

        # Infer the locus trees
        iqtree2 -S concat.best_scheme.nex --prefix loci -T $CORES

        # Compute gene concordance factors
        iqtree2 -te concat.treefile -p concat.best_scheme.nex --gcf loci.treefile --df-tree --cf-verbose --prefix concord -T $CORES

        # Compute site concordance factor using likelihood
        iqtree2 -te concat.treefile -p concat.best_scheme.nex --scfl 100 --prefix concord2 -T $CORES

        # Time divergence estimate
        #   --undo is necessary to continue a previous run when changing/adding options
        echo "Morchella_importuna_SCYDJ1-A1_NA,Verpa_conica_No._21110_NA    -129.6" > calibration_odonnell.txt
        iqtree2 -te concat.treefile -p concat.best_scheme.nex --date calibration_odonnell.txt --date-ci 100 --date-tip 0 --clock-sd 0.2 --undo --prefix divergence_odonnell -T $CORES 
        echo "Choiromyces_venosus_120613-1_NA,Tuber_magnatum_NA_NA    -98.98" > calibration_kraisitudomsook.txt
        echo "Morchella_importuna_SCYDJ1-A1_NA,Verpa_conica_No._21110_NA    -73.2" >> calibration_kraisitudomsook.txt
        iqtree2 -te concat.treefile -p concat.best_scheme.nex --date calibration_kraisitudomsook.txt --date-ci 100 --date-tip 0 --clock-sd 0.2 --undo --prefix divergence_kraisi -T $CORES 
        echo "Choiromyces_venosus_120613-1_NA,Tuber_magnatum_NA_NA    -46.27" > calibration_diaz-escandon.txt
        iqtree2 -te concat.treefile -p concat.best_scheme.nex --date calibration_diaz-escandon.txt --date-ci 100 --date-tip 0 --clock-sd 0.2 --undo --prefix divergence_diaz -T $CORES
    done
done