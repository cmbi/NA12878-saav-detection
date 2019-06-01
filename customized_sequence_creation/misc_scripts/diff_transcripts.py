#!/usr/bin/env python2

import re, collections, sys
from sets import Set

'''This script will help determine how many sequences were affected by the customized sequences
'''

def compare_transcripts(fasta_custom, fasta_reference):
    '''input transcript fastas from the custom and reference'''
    ref=import_fasta(fasta_reference)
    unmatching=Set()
    with open(fasta_custom) as handle:
        seq=''
        for line in handle:
            if line.startswith('>'):
                if seq!='':
                    if '_h' in tid:
                        unmatching.add(tid.split('_')[0])
                    elif tid in ref:
                        if seq!=ref[tid]:
                            unmatching.add(tid)
                    else:
                        raise KeyError('IDs of reference and custom transcript sets dont match')
                    seq=''
                if '|' in line:
                    line=line.split('|')[0]
                tid=line.strip()[1:]
            else:
                seq+=line.strip()
    return len(unmatching)

def import_fasta(fastafile):
    fasta={}
    with open(fastafile) as handle:
        seq=''
        for line in handle:
            if line.startswith('>'):
                if seq!='':
                    fasta[tid]=seq
                    seq=''
                tid=line.strip()[1:]
            else:
                seq+=line.strip()
    return fasta

print compare_transcripts(*sys.argv[1:])