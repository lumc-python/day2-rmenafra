#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 19:31:39 2018

@author: roberta_menafra
"""
import argparse
import gzip
def get_content_gz(fastafile):
    sequence = []    
    with gzip.open(fastafile,'rt') as fh:
        for line in fh:
            if not line.startswith('>'):
                sequence.append(line.strip())
    return sequence
def get_content(fastafile):
    sequence = []    
    with open(fastafile,'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                sequence.append(line.strip())
    return sequence
    
def substrings(sequences, k):
    """Compute substring of n"""
    unique = {}
    for sequence in sequences:
        for i in range(len(sequence)-(k-1)):
            subseq=sequence[i:i+k:]
            if subseq in unique:
                unique[subseq] += 1
            else:
                unique.update({subseq : 1})
    return unique

def calc_gc_percent(sequence):
 
    at_count, gc_count = 0, 0
    for char in sequence.upper():
        if char in ('A', 'T'):
            at_count += 1
        elif char in ('G', 'C'):
            gc_count += 1
        else:
            raise ValueError(
                "Unexpeced character found: {}. Only "
                "ACTGs are allowed.".format(char))
    try:
        return gc_count * 100.0 / (gc_count + at_count)  
    except ZeroDivisionError:
        return 0.0     
##create parser    
parser = argparse.ArgumentParser()
##give argument to the parser
parser.add_argument('input_file', help="input sequence is a fasta file")
parser.add_argument('--gz', action='store_true')
arguments = parser.parse_args()


fastafile = arguments.input_file 
is_compressed = arguments.gz
if is_compressed:
    fastaseq = get_content_gz(fastafile)
else:
    fastaseq = get_content(fastafile)
kmer_dic=substrings(fastaseq,7)

key_max = max(kmer_dic.keys(), key=(lambda x: kmer_dic[x]))
    
print(key_max, calc_gc_percent(key_max))
