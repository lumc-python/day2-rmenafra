# Finding the most common 7-mer in a FASTA file
# =============================================
#
# Write a script to print out the most common 7-mer and its GC percentage from
# all the sequences in data/records.fa. You are free to reuse your existing
# toolbox.
#
# The example FASTA file was adapted from: Genome Biology DNA60 Bioinformatics
# Challenge.
#
# Hints:
# - FASTA files have two types of lines: header lines starting with a ">"
#   character and sequence lines. We are only concerned with the sequence line.
# - Read the string functions documentation.
# - Read the documentation for built in functions.
#
# Challenges:
# - Find out how to change your script so that it can read from
#   data/challenge.fa.gz without unzipping the file first.
# - Can you change the parser so that there is an option flag to tell the
#   program whether the input file is gzipped or not?
# - Can you change your script so that it works for any N-mers instead
#   of for just 7-mers?

import seq_toolbox
import argparse
import gzip

# Dictionary with 7-mer keys and their appearance number as values.
nmer_dict = {}

# Create our argument parser object.
parser = argparse.ArgumentParser()
# Add the expected file path argument.
parser.add_argument('mode', type=str, choices=['textfile', 'gzipfile'],
                    help="Fasta input file path.")
parser.add_argument('file_path', type=str,
                    help="Fasta input file path.")
# Do the actual argument parsing.
args = parser.parse_args()

# Do the file type check.
if args.mode == 'textfile':
    open_any = open
else:
    open_any = gzip.open

# Populate the nmer_dict.
try:
    f = open_any(args.file_path)
    for line in f:
        # We sanitize the file line.
        seq = line.strip().upper()
        if not seq.startswith('>'):
            for i in range(0, len(seq) - 7 + 1):
                nmer_seq = seq[i:i+7]
                if nmer_seq in nmer_dict:
                    nmer_dict[nmer_seq] += 1
                else:
                    nmer_dict[nmer_seq] = 1
finally:
    f.close()

# The most common 7-mer is the nmer_dict key which has the maximum value.
most_nmer = max(nmer_dict, key=nmer_dict.get)

# Get the GC percentage.
gc = seq_toolbox.calc_gc_percent(most_nmer)

# Printing the result.
print "Most common 7-mer is {} with a GC percentage of {:.2f}%.".\
    format(most_nmer, gc)
