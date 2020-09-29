#!/usr/bin/env python2

import re, collections, sys

'''Identify the proteins or RNA that have potential allele specific expression'''

def print_ase_lists(file,txt_general,txt_snv,fa_general,fa_snv,isProtein):
    '''print a list of the identifiers (txt) and the fasta files (fa) of only the allele specific proteins/transcripts
    isProtein=False if rna, True if protein
    '''
    f=open(txt_general,'w')
    d=open(fa_general,'w')
    if isProtein:
        fs=open(txt_snv,'w')
        ds=open(fa_snv,'w')
    seq1=''
    current=''
    add=False
    with open(file) as handle:
        for line in handle:
            if line.startswith(">"):
                if isProtein:
                    if current=='':
                        current=line.split('|')[0]
                        currentfahead=line
                    else:
                        compare=line.split('|')[0]
                        comparefahead=line
                else:
                    if 'haplotype' in line:
                        add=True
                        f.writelines(line.split(' ')[0])
                        d.writelines(line)
                    else:
                        add=False
            else:
                if isProtein:
                    if seq1!='':
                        seq2=line.strip()
                        if seq1!=seq2 and current==compare:
                            mismatch=0
                            if len(seq1)==len(seq2):
                                for idx,aa in enumerate(seq1): #go through and check every amino acid for differences
                                    if aa!=seq2[idx]:
                                        mismatch+=1
                                if mismatch>1: #larger than 1aa difference
                                    f.writelines(current[1:]+'\n')
                                    d.writelines(currentfahead+seq1+'\n'+comparefahead+seq2+'\n')
                                else: #snv
                                    fs.writelines(current[1:]+'\n')
                                    ds.writelines(currentfahead+seq1+'\n'+comparefahead+seq2+'\n')
                            else: #different length sequences, put in general
                                f.writelines(current[1:]+'\n')
                                d.writelines(currentfahead+seq1+'\n'+comparefahead+seq2+'\n')
                        seq1=''
                        current=''
                    else:
                        seq1=line.strip()
                elif add:
                    d.writelines(line)
    return 'output lists written'


print print_ase_lists(sys.argv[1],sys.argv[2],sys.argv[3])