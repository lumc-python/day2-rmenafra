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


def most_common_sequence(f, nmer):
    """
    Determines the most common N-mer sequence from a list of sequences
    present in a fasta file.

    :param f: file handle
    :param nmer: N-mer value
    :return: most common N-mer sequence and its GC percentage.
    """
    # Dictionary with k-mer keys and their appearance number as values.
    nmer_dict = {}
    for line in f:
        # We sanitize the file line.
        seq = line.strip().upper()
        if not seq.startswith('>'):
            for i in range(0, len(seq) - nmer + 1):
                nmer_seq = seq[i:i+nmer]
                if nmer_seq in nmer_dict:
                    nmer_dict[nmer_seq] += 1
                else:
                    nmer_dict[nmer_seq] = 1
    # The most common k-mer is the nmer_dict key which has the maximum value.
    most_nmer = max(nmer_dict, key=nmer_dict.get)

    # Get the GC percentage.
    gc = seq_toolbox.calc_gc_percent(most_nmer)

    return most_nmer, gc

if __name__ == '__main__':
    # Create our argument parser object.
    parser = argparse.ArgumentParser()
    # Add the expected file path argument.
    parser.add_argument('mode', type=str, choices=['textfile', 'gzipfile'],
                        help="Fasta input file path.")
    parser.add_argument('file_path', type=str,
                        help="Fasta input file path.")
    parser.add_argument('nmer', type=int,
                        help="N-mers number.")
    # Do the actual argument parsing.
    args = parser.parse_args()

    # Do the file type check.
    if args.mode == 'textfile':
        open_any = open
    else:
        open_any = gzip.open

    try:
        f = open_any(args.file_path)
        most_nmer, gc = most_common_sequence(f, args.nmer)
    finally:
        f.close()

    # Printing the result.
    print "Most common {}-mer is {} with a GC percentage of {:.2f}%.".\
        format(args.nmer, most_nmer, gc)
