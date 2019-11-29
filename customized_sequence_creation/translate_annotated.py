#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import gffpandas.gffpandas as gffpd
import argparse
import math
from pyteomics.parser import cleave

'''
this script will take the output of the combine_exons script and make translations if a start codon is provided

input: report and transcript fasta output from combine_exons
output: protein fasta of those sequences that had a start codon
'''
def read_in_fasta(fasta):
    df=ph.read_fasta(fasta)
    df=df[['id','sequence']]
    # if df['id'].str.contains('|').all():
    #     df['id']=df['id'].str.split('|').apply(lambda x: x[0]) #keep only the name of the transcript
    return(df)

def translate(dna):
    protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V", \
    "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V", "TTA" : "L", \
     "CTA" : "L", "ATA" : "I", "GTA" : "V","TTG" : "L", "CTG" : "L", \
     "ATG" : "M", "GTG" : "V","TCT" : "S", "CCT" : "P", "ACT" : "T", \
     "GCT" : "A","TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A", \
     "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A","TCG" : "S", \
     "CCG" : "P", "ACG" : "T", "GCG" : "A","TAT" : "Y", "CAT" : "H", \
     "AAT" : "N", "GAT" : "D", "TAC" : "Y", "CAC" : "H", "AAC" : "N", \
     "GAC" : "D","TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E", \
     "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E", "TGT" : "C", \
     "CGT" : "R", "AGT" : "S", "GGT" : "G","TGC" : "C", "CGC" : "R", \
     "AGC" : "S", "GGC" : "G","TGA" : "STOP", "CGA" : "R", "AGA" : "R", \
     "GGA" : "G","TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"}
    protein_sequence = ""
    # Generate protein sequence
    for i in range(0, len(dna)-(3+len(dna)%3), 3):
        if protein[dna[i:i+3]] == "STOP" :
            break
        protein_sequence += protein[dna[i:i+3]]
    return(protein_sequence)

def translate_from_startcodon(sequence,startpos):
    '''translate sequence from specified start codon position'''
    return(translate(sequence[int(startpos):]))

def renumber_variants_from_startcodon(variant_list,start):
    '''get the protein coordinate of the variants'''
    return([math.ceil((int(x)-int(start))/3) for x in variant_list])

def digest(protein_sequence,het_variants,het_origins,hom_variants,hom_origins):
    '''digests the protein sequences with trypsin, fetches the start position in the sequence and then checks if 
    a variant position falls in between
    returns list of "variant info", list of strings
    each string is a comma seperated list of information
    ['peptide','peptide_start','var_pos','var_org','var_type']
    '''
    seq_cut = cleave(protein_sequence, '[KR]', 2)
    variant_info=[] #the list to return
    plist=set() #to prevent duplicates
    for peptide in seq_cut:
        if peptide in plist:
            continue
        plist.add(peptide)
        if len(peptide) < 6: #only allow peptides of at least 6 aa
            continue
        start = find_str(protein_sequence, peptide)
        # for variant_list in [het_variants,hom_variants]
        for i,l in enumerate(het_variants):
            if start < int(l) < start+len(peptide):
                variant_info.append(',',join([peptide,start,l,het_origins[i],'heterozygous']))
        for ind,var in enumerate(hom_variants):
            if start < int(var) < start+len(peptide):
                variant_info.append(',',join([peptide,start,var,hom_origins[ind],'homozygous']))
    return(variant_info)


def find_str(s, char):
    index = 0
    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index+len(char)] == char:
                    return index
            index += 1
    return -1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--fasta', help='Transcript fasta', required=True)
    parser.add_argument('--report', help='File report', required=True)
    args=vars(parser.parse_args())
    #import the fasta
    seqs=read_in_fasta(args['fasta'])
    #import the information
    var_info=pd.read_csv(args['report'],sep='\t',header=0,names=['id','pos_het','het_origin','pos_hom','hom_origin','start_codon_position'])
    #only include those with a start codon
    var_info=var_info[var_info['start_codon_position']!='None']
    #merge
    df=pd.merge(var_info,seqs,on='id')
    #unstack multiple start codons- so only one start codon position per row
    df['start_codon_position']=df['start_codon_position'].str.split(',') #split start codon positions in case there is more than one
    unstacked=df.apply(lambda x: pd.Series(x['start_codon_position']),axis=1).stack().reset_index(level=1, drop=True))
    unstacked.name='start_codon_position'
    df=df.drop('start_codon_position').join(unstacked)
    #renumber the variant positions
    for colname in ['pos_het','het_origin','pos_hom','hom_origin']:
        df[colname]=df[colname].str.split(',')
        if colname in ['pos_het','pos_hom']:
            df[colname]=df.apply(lambda x: renumber_variants_from_startcodon(x[colname],x['start_codon_position']),axis=1)
    #predict protein sequence
    df['protein_sequence']=df.apply(lambda x: translate_from_startcodon(x['sequence'],x['start_codon_position']),axis=1)
    #digest protein sequence and get peptides with variants in them
    df['pep_to_variant']=df.apply(lambda x:)
