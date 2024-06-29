#!/usr/bin/env python
import os
import argparse
import sys
import glob
from Bio import SeqIO

# Version
_verion_= "0.1"

# Argparse Setup
parser = argparse.ArgumentParser(description="A tool for extracting genes from a genbank file and writing to individual fastas.")
parser.add_argument("-i", "--input", required=True, help="Path to file containing locus tags of interest.")
parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file.")
parser.add_argument("-o", "--output", required=True, help="Path to output directory.")
args = parser.parse_args()

# Welcome
print('###################')
print('Welcome to GBKsplit')
print('###################')
print('')

# Directory orientation
invoked_from = os.getcwd()
if not os.path.exists(args.output):
    os.mkdir(args.output)
    print('Output directory created.','\n')
os.chdir(args.output)
output_directory = os.getcwd()
os.chdir(invoked_from)

# Loads locus tags from file and creates output filenames.
try:
    locus_tags = open(args.input).read().splitlines()
    if len(locus_tags) == 0:
        print('Could not find any locus tags in input file.')
        print('')
        print('Please check your file and try again.')
        print('')
        sys.exit(1)
    print('Loaded', len(locus_tags),'locus tags.','\n','\n')
except Exception as e:
    print('Error reading input file.', '\n')
    print(e,'\n')
    sys.exit(1)
output_names = [x + '.fasta' for x in locus_tags]
output_files = [output_directory + '/' + x for x in output_names]

# Process genbank file:
print('Commencing processing of genbank file.','\n')
if len(locus_tags) >= 100:
    print('This make take some time...','\n')
for input_tags, output_fasta in zip(locus_tags, output_files):
    with open(output_fasta, 'w') as out:
        for entries in SeqIO.parse(args.genbank,"genbank"):
            if entries.features:
                for feature in entries.features:
                    if feature.type == "CDS":
                        if feature.qualifiers['locus_tag'][0] == input_tags:
                            out.write(">%s\n%s\n" % (
                                feature.qualifiers['locus_tag'][0],
                                feature.location.extract(entries).seq))

# Farewell
print('Done!','\n')
print('Author: www.github.com/stevenjdunn','\n','\n')
print('#########')
print('Finished!')
print('#########')
