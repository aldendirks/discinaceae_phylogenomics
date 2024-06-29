# _Discinaceae_ phylogenomics pipeline

<!-- vscode-markdown-toc -->
* 1. [Setup](#Setup)
	* 1.1. [Software](#Software)
	* 1.2. [Set up directory structure](#Setupdirectorystructure)
* 2. [Genome assembly and phylogenomics](#Genomeassemblyandphylogenomics)
	* 2.1. [Download Illumina data from NCBI SRA](#DownloadIlluminadatafromNCBISRA)
	* 2.2. [Check read quality and trim adapters](#Checkreadqualityandtrimadapters)
	* 2.3. [Correct reads and assemble](#Correctreadsandassemble)
	* 2.4. [Filter ascomata genomes](#Filterascomatagenomes)
	* 2.5. [Rescue *Gyromitra leucoxantha* MICH25407](#RescueGyromitraleucoxanthaMICH25407)
	* 2.6. [Further contaminant filtering with NCBI Foreign Contamination Screen](#FurthercontaminantfilteringwithNCBIForeignContaminationScreen)
	* 2.7. [TO DO: Organize genomes](#TODO:Organizegenomes)
	* 2.8. [TO DO: Download outgroup genomes](#TODO:Downloadoutgroupgenomes)
	* 2.9. [Calculate genome coverage](#Calculategenomecoverage)
	* 2.10. [Extract rDNA sequences](#ExtractrDNAsequences)
	* 2.11. [Run BUSCO](#RunBUSCO)
	* 2.12. [Make phylogenomic tree](#Makephylogenomictree)
* 3. [ Genome annotation and summary statistics](#Genomeannotationandsummarystatistics)
	* 3.1. [Predict genes](#Predictgenes)
	* 3.2. [Identify secondary metabolite biosynthesis gene clusters](#Identifysecondarymetabolitebiosynthesisgeneclusters)
	* 3.3. [Assign gene function](#Assigngenefunction)
	* 3.4. [Complete annotation](#Completeannotation)
	* 3.5. [Compile genome stats](#Compilegenomestats)
* 4. [Mating system analysis](#Matingsystemanalysis)
	* 4.1. [Soft link proteomes into one directory](#Softlinkproteomesintoonedirectory)
	* 4.2. [Make proteomes BLAST database](#MakeproteomesBLASTdatabase)
	* 4.3. [Identify MAT genes and create clinker map](#IdentifyMATgenesandcreateclinkermap)
* 5. [TO DO: Figures and statistical analysis](#TODO:Figuresandstatisticalanalysis)
	* 5.1. [TO DO: Copy data](#TODO:Copydata)
* 6. [Supplementary](#Supplementary)
	* 6.1. [Transfer data from UM Advanced Genomics Core to Great Lakes computing cluster](#TransferdatafromUMAdvancedGenomicsCoretoGreatLakescomputingcluster)
		* 6.1.1. [Globus](#Globus)
		* 6.1.2. [Check the demultiplexing statistics (Illumina)](#CheckthedemultiplexingstatisticsIllumina)
		* 6.1.3. [Check file completeness](#Checkfilecompleteness)
	* 6.2. [Upload data to NCBI](#UploaddatatoNCBI)
		* 6.2.1. [BioProject](#BioProject)
		* 6.2.2. [BioSample](#BioSample)
		* 6.2.3. [SRA](#SRA)
		* 6.2.4. [Genomes](#Genomes)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

##  1. <a name='Setup'></a>Setup

###  1.1. <a name='Software'></a>Software

The following software need to be installed and findable in $PATH.

**Published software used in this pipeline**
* `antiSMASH` 7.1.0 and dependencies (including meme 4.11.2 for RiPP detection) contained in a `Conda` environment called `antismash`
* `ASTRAL` 5.15.5
* `BBTools` 38.96
* `bedtools` 2.30.0
* `BUSCO` 5.5.0
* `Clinker` 0.0.28
* `CUDA` 11.5.1 (for computing on GPUs with SCGid)
* `FastQC` 0.11.9
* `FastTree` 2.1.10
* `Funannotate` 1.8.13 and dependencies contained in a `Conda` environment called `funannotate`
* `IQTREE2` 2.2.6 (built with LSD2 for time divergence estimation)
* `ITSx` 1.1.3
* `MUSCLE` 3.8.31
* `NCBI BLAST+` 2.12.0
* `Open MPI` 4.1.4 (for parallel processing with SCGid)
* `Python` 3.10.4
* `QUAST` 5.2.0
* `R` 4.2.0
* `SCGid` 0.9b0
* `Seqtk` 1.3
* `SPAdes` 3.15.5
* `SRA Toolkit` 2.10.9
* `trimAl` 1.4.rev22
* `Trimmomatic` 0.36

**Other software used in this pipeline**
* [BUSCO phylogenomics pipeline](https://github.com/jamiemcg/BUSCO_phylogenomics/releases/tag/3), BUSCO V3 (Old) release - have `BUSCO_phylogenomics.py` findable in `PATH`
* [Extract-ITS-sequences-from-a-fungal-genome](https://github.com/aldendirks/Extract-ITS-sequences-from-a-fungal-genome), patch-1 branch, forked from pwkooij/Extract-ITS-sequences-from-a-fungal-genome, forked from fantin-mesny/Extract-ITS-sequences-from-a-fungal-genome - have `extractITS.py` findable in `PATH`

###  1.2. <a name='Setupdirectorystructure'></a>Set up directory structure

The directory containing this README.md file will be the project directory. Change to the parent directory and run the code below to clone the GitHub directory `discinaceae_phylogenomics`. 

```
git clone https://github.com/aldendirks/discinaceae_phylogenomics.git
export PROJECT_DIR="$PWD/discinaceae_phylogenomics"
cd $PROJECT_DIR
mkdir -p antismash busco coverage blastdb data fastqc funannotate genomes itsx mat phylogenomics quast rdna scgid seqs spades summary
```

Make sure the list of samples is unix compatible and ends with a space.

```
export SAMPLE_LIST="$PROJECT_DIR/misc/samples.tsv"
dos2unix $SAMPLE_LIST
function file_ends_with_newline() {
    [[ $(tail -c1 "$1" | wc -l) -gt 0 ]]
}
if ! file_ends_with_newline $SAMPLE_LIST
then
    echo "" >> $SAMPLE_LIST
fi
```

##  2. <a name='Genomeassemblyandphylogenomics'></a>Genome assembly and phylogenomics

###  2.1. <a name='DownloadIlluminadatafromNCBISRA'></a>Download Illumina data from NCBI SRA

Skip to step 2.7 if you'd like to go right to downloading the assembled genomes from NCBI and preparing them for phylogenomics analysis.  

Download raw reads for processing. Further information on downloading data from NCBI SRA using `SRA Toolkit` can be found [here](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump).

```
bash $PROJECT_DIR/scripts/00_download_sra_data.sh
```

###  2.2. <a name='Checkreadqualityandtrimadapters'></a>Check read quality and trim adapters

```
sbatch $PROJECT_DIR/scripts/01_check_read_quality_and_trim.sh
```

###  2.3. <a name='Correctreadsandassemble'></a>Correct reads and assemble

First, run read error correction for ascomata samples. This can be a memory-intensive step and might require adjustment for different samples (memory requirement up to 200 GB). By doing this first, SPAdes assembly can proceed in a more predictable fashion. 

```
while IFS=$'\t' read -r -a samples_array
do
    if [ "${samples_array[7]}" = "ingroup" ] && [ "${samples_array[14]}" != "culture" ]; then
        BATCH="${samples_array[1]}"
        SAMPLE="${samples_array[2]}"
        sbatch $PROJECT_DIR/scripts/02_correct_read_errors.sh $BATCH $SAMPLE
    fi
done < $SAMPLE_LIST
```

Then assemble genomes with `SPAdes`.

```
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ]; then
        if [ "${samples_array[14]}" = "culture" ]; then
            sbatch $PROJECT_DIR/scripts/03a_assemble_genome_isolate.sh $BATCH $SAMPLE
        else
            sbatch $PROJECT_DIR/scripts/03b_assemble_genome_meta.sh $BATCH $SAMPLE
        fi
    fi
done < $SAMPLE_LIST
```

###  2.4. <a name='Filterascomatagenomes'></a>Filter ascomata genomes

SCGid removes contaminant contigs by triangulating evidence from three different modules. Make sure you have activated the `SCGid` virtual environment before running the jobs with slurm. 

First, run the GCT module. 

```
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ] && [ "${samples_array[14]}" != "culture" ]; then
        sbatch $PROJECT_DIR/scripts/04a_run_scgid_gct.sh $BATCH $SAMPLE
    fi
done < $SAMPLE_LIST
```

Then, the codons module. 

```
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ] && [ "${samples_array[14]}" != "culture" ]; then
        sbatch $PROJECT_DIR/scripts/04b_run_scgid_codons.sh $BATCH $SAMPLE
    fi
done < $SAMPLE_LIST
```

Next, the first part of the kmers module.

```
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ] && [ "${samples_array[14]}" != "culture" ]; then
        sbatch $PROJECT_DIR/scripts/04c_run_scgid_kmers.sh $BATCH $SAMPLE
    fi
done < $SAMPLE_LIST
```

Manually select regions from the ESOM map. See the [`SCGid` GitHub page](https://github.com/amsesk/SCGid) for more information. Some additions to the notes there:

* A large `.wts` file may take a while to load. 
* After loading the `.wts` file, you also need to load the `.bm` file. 
* It is useful to select "Tiled display" in the ESOM window and change Zoom to 1x. 
* Under the "Classes" tab, I like to change the color of Fungi to something that pops and is easy to trace, like neon yellow. Change the color of any class by selecting the color next to the text. 
* Bold the Fungi dots so they are easier to see by clicking the box on the left. 
* Right click to finalize your the selection after you click around the fungi blobs. 
* I'd recommend using an external mouse to make your selections. It is easy to make an accidental right click when you are clicking with a trackpad. 

Finally, finish the kmers module and generate a consensus sequence. Note, the ESOM selected class ID is set to "6" within the script. While this should usually be the correct class, you might need to change it so it matches your selection. 

```
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ] && [ "${samples_array[14]}" != "culture" ]; then
        sbatch $PROJECT_DIR/scripts/04d_run_scgid_kmers_extract.sh $BATCH $SAMPLE
    fi
done < $SAMPLE_LIST
```

###  2.5. <a name='RescueGyromitraleucoxanthaMICH25407'></a>Rescue *Gyromitra leucoxantha* MICH25407 

This genome is the only one with significant fungal contamination. `SCGid` does not filter out the contaminant fungal genome. The *Gyromitra leucoxantha* genome is extracted by BLASTing the raw reads against the assembled and filtered *G. leucoxantha* MICH352087 genome and using those hits in assembly. 

```
sbatch $PROJECT_DIR/scripts/05_rescue_gyromitra_leucoxantha.sh
```

###  2.6. <a name='FurthercontaminantfilteringwithNCBIForeignContaminationScreen'></a>Further contaminant filtering with NCBI Foreign Contamination Screen

Upon submission to NCBI, it was found that `SCGid` did not completely filter out contaminant contigs according to NCBI's Foreign Contamination Screen (FCS) tool, which is automatically run on submitted genomes. The automatically filtered genomes were downloaded and used in subsequent steps. The FCS tool suite can be installed locally but requires very high memory and large harddrive storage to run efficiently. It may be that filtering programs like `SCGid` can be skipped entirely and contaminat filtering done automatically by uploading genomes to NCBI. Genomes need to pass this screening step to be uploaded to NCBI anyway so one might as well us their tool. 

###  2.7. <a name='TODO:Organizegenomes'></a>TO DO: Organize genomes

TO DO: Get GenBank accession numbers and make sure they can be downloaded

Skip these first two blocks of code if you would like to go right to downloading the assembled genomes. 

Softlink the assembled genomes into the `genomes` directory with a new name.

```
export GENOMES_DIR="$PROJECT_DIR/genomes"
while IFS=$'\t' read -r -a samples_array
do
    BATCH="${samples_array[1]}"
    SAMPLE="${samples_array[2]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')" # substitute spaces with underscores and remove any trailing spaces or underscores
    if [ "$BATCH" != "JGI" ] && [ "$BATCH" != "NCBI" ]; then
        if [ "${samples_array[14]}" = "culture" ]; then
            FILE="$PROJECT_DIR/spades/Sample_${BATCH}-AD-${SAMPLE}-isolate/scaffolds_trimmed.fasta"
        else
            FILE="$PROJECT_DIR/scgid/Sample_${BATCH}-AD-${SAMPLE}-meta/meta_scgid_output/consensus/meta.consensus.filtered.assembly.fasta"
        fi
    fi
    ln -s $FILE $GENOMES_DIR/$NAME.fasta
done < $SAMPLE_LIST
```

Download the JGI MycoCosm *Gyromitra* genomes. You will need to have a profile set up with JGI. Replace "[email_address]" and "[password]" (brackets inclusive) with your JGI login email and password, respectively. 

```
# Log in to JGI
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=[email_address]' --data-urlencode 'password=[password]' -c ~/cookies > /dev/null

# Download genomes 
while IFS=$'\t' read -r -a samples_array
do
    if [ "${samples_array[7]}" = "ingroup" ] && [ "${samples_array[1]}" = "JGI" ]; then
        PORTAL="${samples_array[2]}"
        NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
        curl "https://genome.jgi.doe.gov/portal/$PORTAL/download/${PORTAL}_AssemblyScaffolds.fasta.gz" -b ~/cookies | gunzip > $GENOMES_DIR/$NAME.fasta
    fi
done < $SAMPLE_LIST
```

Conversely, if you skipped genome assembly, you can download the assembled genomes from JGI and NCBI. See the instructions for the previous block of code on how to download genomes from JGI.  

```
export GENOMES_DIR="$PROJECT_DIR/genomes"
while IFS=$'\t' read -r -a samples_array
do
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "${samples_array[7]}" = "ingroup" ]; then
        if [ "${samples_array[1]}" = "JGI" ]; then
            PORTAL="${samples_array[2]}"
            curl "https://genome.jgi.doe.gov/portal/$PORTAL/download/${PORTAL}_AssemblyScaffolds.fasta.gz" -b ~/cookies | gunzip > $GENOMES_DIR/$NAME.fasta
        else
            ACCESSION="${samples_array[6]}"
            curl "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$ACCESSION/download?include_annotation_type=GENOME_FASTA" > $GENOMES_DIR/$NAME
            unzip $GENOMES_DIR/$NAME
            mv $GENOMES_DIR/ncbi_dataset/data/$ACCESSION/*.fna $GENOMES_DIR/$NAME.fasta
            rm -r $GENOMES_DIR/ncbi_dataset $GENOMES_DIR/$NAME $GENOMES_DIR/README.md
        fi
    fi
done < $SAMPLE_LIST
```

###  2.8. <a name='TODO:Downloadoutgroupgenomes'></a>TO DO: Download outgroup genomes

TO DO: Maybe change 7th column (WGS accession) to GenBank accessions when I have it and change 2 for samples_array NCBI to 6. 

```
while IFS=$'\t' read -r -a samples_array
do
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "${samples_array[7]}" = "outgroup" ]; then
        if [ "${samples_array[1]}" = "JGI" ]; then
            PORTAL="${samples_array[2]}"
            curl "https://genome.jgi.doe.gov/portal/$PORTAL/download/${PORTAL}_AssemblyScaffolds.fasta.gz" -b ~/cookies | gunzip > $GENOMES_DIR/$NAME.fasta
        elif [ "${samples_array[1]}" = "NCBI" ]; then
            ACCESSION="${samples_array[2]}"
            curl "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$ACCESSION/download?include_annotation_type=GENOME_FASTA" > $GENOMES_DIR/$NAME
            unzip $GENOMES_DIR/$NAME
            mv $GENOMES_DIR/ncbi_dataset/data/$ACCESSION/*.fna $GENOMES_DIR/$NAME.fasta
            rm -r $GENOMES_DIR/ncbi_dataset $GENOMES_DIR/$NAME $GENOMES_DIR/README.md
        fi
    fi
done < $SAMPLE_LIST
```

###  2.9. <a name='Calculategenomecoverage'></a>Calculate genome coverage

```
sbatch $PROJECT_DIR/scripts/09_calculate_coverage.sh
```

###  2.10. <a name='ExtractrDNAsequences'></a>Extract rDNA sequences

```
sbatch $PROJECT_DIR/scripts/10_extract_rdna.sh
```

###  2.11. <a name='RunBUSCO'></a>Run BUSCO

BUSCO stands for benchmarking universal single-copy orthologs. BUSCOs are used in the phylogenomic analysis. 

```
cd $PROJECT_DIR/busco
busco --download ascomycota_odb10
for SAMPLE in $PROJECT_DIR/genomes/*
do
    sbatch $PROJECT_DIR/scripts/06_run_busco.sh $SAMPLE
done
```

###  2.12. <a name='Makephylogenomictree'></a>Make phylogenomic tree

The BUSCO phylogenomics pipeline is used for coalescence and to prepare the dataset for `IQTREE2` concatenation analysis. You will need to define the path to the `ASTRAL` library files and `ASTRAL` jar file, e.g., `"/home/adirks/apps/ASTRAL/Astral/lib/"` and `"/home/adirks/apps/ASTRAL/Astral/astral.5.15.5.jar"`.

```
export ASTRAL_LIB_PATH="/path/to/Astral/lib/"
export ASTRAL_JAR_PATH="/path/to/astral.5.15.5.jar"
# consider running the following code to find the paths if available, although find can take a while
# find / -name "astral.5.15.5.jar"

sbatch $PROJECT_DIR/scripts/07_phylogenomics.sh
```

##  3. <a name='Genomeannotationandsummarystatistics'></a> Genome annotation and summary statistics

Genome annotation is completed with `Funannotate`, available through [GitHub](https://github.com/nextgenusfs/funannotate) and documented [here](https://funannotate.readthedocs.io/en/latest/). 

###  3.1. <a name='Predictgenes'></a>Predict genes

First, predict genes. The JGI <i>Gyromitra esculenta</i> expressed sequence tags (EST) file is used to help with gene prediction. 

```
conda deactivate
conda activate funannotate

# Get EST file from JGI
curl "https://genome.jgi.doe.gov/portal/Gyresc1/download/Gyresc1_ESTs_20160326_est.fasta.gz" -b ~/cookies | gunzip > $PROJECT_DIR/funannotate/Gyresc1_ESTs_20160326_est.fasta

# Run Funannotate
while IFS=$'\t' read -r -a samples_array
do
    SPECIES="${samples_array[9]}"
    ACCESSION="${samples_array[11]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "${samples_array[7]}" = "ingroup" ]; then
        sbatch $PROJECT_DIR/scripts/08a_run_funannotate_predict.sh -s "'"$SPECIES"'" -i $ACCESSION -g $GENOMES_DIR/$NAME.fasta
    fi
done < $SAMPLE_LIST
```

At this point, an NCBI bioproject number, individual biosample numbers, and a submission template (`.sbt` file) were generated before continuing with annotation. The `template.sbt` file I used is available in the `misc` directory.

###  3.2. <a name='Identifysecondarymetabolitebiosynthesisgeneclusters'></a>Identify secondary metabolite biosynthesis gene clusters

Identify secondary metabolism clusters with `antiSMASH`.

```
conda deactivate
conda activate antismash
while IFS=$'\t' read -r -a samples_array
do
    SPECIES="${samples_array[9]}"
    ACCESSION="${samples_array[11]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "${samples_array[7]}" = "ingroup" ]; then
        sbatch $PROJECT_DIR/scripts/08b_run_antismash7.sh -s "'"$SPECIES"'" -i $ACCESSION -g $GENOMES_DIR/$NAME.fasta
    fi
done < $SAMPLE_LIST
```

###  3.3. <a name='Assigngenefunction'></a>Assign gene function

Use `InterProScan` for functional analysis of proteins.

```
conda deactivate
conda activate funannotate
while IFS=$'\t' read -r -a samples_array
do
    SPECIES="${samples_array[9]}"
    ACCESSION="${samples_array[11]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "${samples_array[7]}" = "ingroup" ]; then
        sbatch $PROJECT_DIR/scripts/08c_run_interproscan.sh -s "'"$SPECIES"'" -i $ACCESSION -g $GENOMES_DIR/$NAME.fasta
    fi
done < $SAMPLE_LIST
```

###  3.4. <a name='Completeannotation'></a>Complete annotation

```
while IFS=$'\t' read -r -a samples_array
do
    SPECIES="${samples_array[9]}"
    ACCESSION="${samples_array[11]}"
    LOCUS_TAG="${samples_array[3]}"
    NAME="$(echo ${samples_array[15]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    if [ "${samples_array[7]}" = "ingroup" ]; then
        sbatch $PROJECT_DIR/scripts/08d_run_funannotate_annotate.sh -g $GENOMES_DIR/$NAME.fasta -l $LOCUS_TAG
    fi
done < $SAMPLE_LIST
```

###  3.5. <a name='Compilegenomestats'></a>Compile genome stats

Get genome stats.

```
bash $PROJECT_DIR/scripts/11a_summarize_genome_stats.sh
```

Get `antiSMASH` secondary metabolite gene cluster counts.

```
bash $PROJECT_DIR/scripts/11b_summarize_antismash_counts.sh
```

Get CAZyme family counts. 

```
bash $PROJECT_DIR/scripts/11c_summarize_cazyme_counts.sh
```

##  4. <a name='Matingsystemanalysis'></a>Mating system analysis 

###  4.1. <a name='Softlinkproteomesintoonedirectory'></a>Soft link proteomes into one directory

```
bash $PROJECT_DIR/12_ln_proteomes.sh
```

###  4.2. <a name='MakeproteomesBLASTdatabase'></a>Make proteomes BLAST database

```
bash $PROJECT_DIR/13_make_blast_db.sh
```

###  4.3. <a name='IdentifyMATgenesandcreateclinkermap'></a>Identify MAT genes and create clinker map

```
bash $PROJECT_DIR/14_blast_mat_genes.sh
```

##  5. <a name='TODO:Figuresandstatisticalanalysis'></a>TO DO: Figures and statistical analysis

###  5.1. <a name='TODO:Copydata'></a>TO DO: Copy data

Get data into one directory.

```
cp $PROJECT_DIR/summary/counts_antismash.txt $PROJECT_DIR/data
cp $PROJECT_DIR/summary/counts_cazy.txt $PROJECT_DIR/data
cp $PROJECT_DIR/summary/genome_stats.txt $PROJECT_DIR/data
```

##  6. <a name='Supplementary'></a>Supplementary

###  6.1. <a name='TransferdatafromUMAdvancedGenomicsCoretoGreatLakescomputingcluster'></a>Transfer data from UM Advanced Genomics Core to Great Lakes computing cluster

####  6.1.1. <a name='Globus'></a>Globus

Click on the Globus link in the email you receive from the Advanced Genomics Core data team with your data. This will take you to Globus File Manager. You should see two panels. The left panel should be the Core's Globus endpoint and should show your data.

On the right panel, type or select umich#greatlakes as the "Collection". This should then take you to your home directory in Great Lakes. Enter "/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/" as the "Path". Create or navigate to a directory to store your raw data.

On the left panel, select all of your data. Click "Start" to begin the data transfer. You can close the window; the transfer happens in the background. You can check the status of the transfer under the "Activity" button on the far left menu in Globus. A dataset of around 50 GB should only take a few minutes to complete. 

Please add a plain text README file in the directory that describes what your data are.

####  6.1.2. <a name='CheckthedemultiplexingstatisticsIllumina'></a>Check the demultiplexing statistics (Illumina)

The file `DemuxStats.csv` is included with Illumina data and shows your samples various information, including the number of reads. You multiple your reads by the read length to get your total number of base pairs per sample. If you have an estimate for the total genome size, you can get a sense if your coverage goals were met. 

####  6.1.3. <a name='Checkfilecompleteness'></a>Check file completeness

This website has information on checking file completeness using the md5 file: https://brcf.medicine.umich.edu/cores/advanced-genomics/data-delivery/

Under Frequently Asked Questions, "How do I know I have a complete copy of the data?", it says:

```
cd <service-request-name>
md5sum –-quiet --check <service-request-name>.md5
```

The `--quiet` option means it will only show output if there is a problem, so the output will be blank and processing for a while if everything is OK.

###  6.2. <a name='UploaddatatoNCBI'></a>Upload data to NCBI

####  6.2.1. <a name='BioProject'></a>BioProject

Create a [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/) to receive an accession number that will be used to link all your data. The BioProject submission portal is [here](https://submit.ncbi.nlm.nih.gov/subs/bioproject/).

> A BioProject is a collection of biological data related to a single initiative, originating from a single organization or from a consortium. A BioProject record provides users a single place to find links to the diverse data types generated for that project.

####  6.2.2. <a name='BioSample'></a>BioSample

Once you have a BioProject accession, proceed to [creating BioSamples](https://submit.ncbi.nlm.nih.gov/subs/biosample/). For a fungal genome, use the "MIGS Eukaryotic" data template package, even if the genome was filtered from a mushroom and thus the raw data are more like a metagenome. Locus tag prefixes are generated for each BioSample. You should receive an email within a few days with those numbers, which then become available on the BioProject page. `Funannotate` uses locus tags as prefixes in the FASTA headers, which are required for annotated genome submission. 

> The BioSample database contains descriptions of biological source materials used in experimental assays.

####  6.2.3. <a name='SRA'></a>SRA

Once you have BioSample accessions, proceed to the [NCBI SRA submission portal](https://submit.ncbi.nlm.nih.gov/subs/sra/) to get SRA accession numbers and to upload the raw reads. Do not modify the name of the SRA template file that you download from NCBI.

> Sequence Read Archive (SRA) data, available through multiple cloud providers and NCBI servers, is the largest publicly available repository of high throughput sequencing data. The archive accepts data from all branches of life as well as metagenomic and environmental surveys. SRA stores raw sequencing data and alignment information to enhance reproducibility and facilitate new discoveries through data analysis.

I recommend uploading the raw read data to NCBI from Great Lakes via FTP. Using this option, you will need to put all the files in a single directory, for example one named `ncbi_sra` in your scratch. The file names need to match what you included in the SRA data file. 

Follow these instructions to upload your data:
* Navigate to the local directory with your data
* Establish an FTP connection with `ftp -i`
* Type `open ftp-private.ncbi.nlm.nih.gov`
* Provide the username listed in the Submission Portal, likely `subftp`
* Enter the password listed in the Submission Portal (unique to you)
* Navigate to your account folder listed in the Sumbission Portal (unique to you)
* Make a new directory (e.g., `mkdir sra_submission`)
* Navigate to that direcotry (`cd sra_submission`)
* Transfer all your local files to this directory with `mput *`
* Go back to the Submission Portal page and select the upload folder. It takes at least 10 minutes for transferred files to appear in the preload option.
* Exit FTP with the command `bye`

####  6.2.4. <a name='Genomes'></a>Genomes

Submit genomes [here](https://submit.ncbi.nlm.nih.gov/about/genome/). You are going to upload a `.fasta` file for an unannotated genome or a `.sqn` file for an annotated genome (generated by `Funannotate`). Follow the instructions above for SRA files to upload the genome files. Make sure to put all your `.fasta` or `.sqn` files together in the same directory. I recommend only uploading the genomes (unannotated FASTA file) because NCBI will probably identify some issues with the annotation that might be difficult to resolve. 