# GBKSplit
**Extracts genes using a list of locus tags from .gbk to logically named nucleotide .fasta's.**

## Dependencies
- Python 3
- Biopython

## Quick start
    gbksplit.py -i /path/to/input/locus_tags.txt -g /path/to/genbank/file.gbk -o /path/to/output_directory/
    
## Usage:     
    usage: gbksplit.py [-h] -i INPUT -g GENBANK -o OUTPUT

    A tool for splitting up .fasta files, creates a new file for each contig.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Path to file containing locus tags of interest.
      -g GENBANK, --genbank GENBANK
                            Path to genbank file.
      -o OUTPUT, --output OUTPUT
                            Path to output directory.
## How
Simply point gbksplit to a text file containing a list of locus tags, separated by a new line. The script uses Biopython's SeqIO to match this to the CDS record and export the nucleotide sequence in separate .fasta files named after the provided locus_tag. 

## Why?
I use it as part of a SNP calling pipeline to extract polymorphic genes for further downstream analysis. 

## Notes
Takes 1-3 seconds per locus tag on a 1.2 Ghz core m3, depending on (I think) gene size.
