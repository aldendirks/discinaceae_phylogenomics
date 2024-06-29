"""
Modified from: https://bioinformatics.stackexchange.com/questions/4365/how-to-extract-the-protein-fasta-file-from-a-genbank-file
Usage: gbk_converter.py foo.gbk
"""

import sys, os
from Bio import SeqIO

file_name = sys.argv[1]
header_name = sys.argv[2]

# stores all the CDS entries
all_entries = []

with open(file_name, 'r') as GBFile:

    GBcds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)

    for cds in GBcds:
        if cds.seq is not None:
            # print(cds) # print record to see all the availale attributes and annotations
            cds.id = header_name
            cds.description = cds.annotations["gene"]
            all_entries.append(cds)


# write file
SeqIO.write(all_entries, '{}.fasta'.format(file_name[:-4]), 'fasta')