#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import gffpandas.gffpandas as gffpd
import argparse
import math
import phylopandas as ph
from pyteomics.parser import cleave

'''
this script will take the output of the combine_exons script and make translations
also spits out all peptides that correspond with a variant position (whether or not synonymous)

input: report and transcript fasta output from combine_exons, angel .cds if novel sequences
output: protein fasta, csv file of peptides containing a variant position
'''
def read_in_fasta(fasta,cds=False):
    '''read fasta'''
    df=ph.read_fasta(fasta)
    df['id']=df['id'].apply(separate_haplotype)
    df[['id','haplotype']]=pd.DataFrame(df['id'].values.tolist(),index=df.index)
    if cds:
        df['start_codon_position']=df['description'].str.split('pos:').apply(lambda x: x[1])
        df['start_codon_position']=df['start_codon_position'].str.split('-').apply(lambda x: x[0])
        return(df[['id','haplotype','start_codon_position']])
    return(df[['id','haplotype','sequence']])

def write_to_fasta(df,outputfile):
    assert len(df.columns)==3, f"wrong number of columns: {str(len(df.columns))}"
    df.columns=['id','haplotype', 'sequence']
    df['id']='>'+df['id']+' hap:'+df['haplotype']
    df[['id','sequence']].to_csv(outputfile,header=None,index=None,sep='\n')

def separate_haplotype(tid):
    if '|' in tid:
        tid=tid.split('|')[0]
    if '_h' in tid:
        return(tid.split('_h'))
    return([tid,'None'])

def translate(dna):
    ''''''
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

def renumber_variants(variant_list):
    '''get the protein coordinate of the variants'''
    if 'None' not in variant_list[0]:
        return([math.ceil(int(x)/3) for x in variant_list])
    return(variant_list)

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
            if isinstance(l,int) and start <= l and l <= start+len(peptide):
                variant_info.append(','.join([peptide,str(start),str(l),het_origins[i],'heterozygous']))
        for ind,var in enumerate(hom_variants):
            if isinstance(var,int) and start <= var and var <= start+len(peptide):
                variant_info.append(','.join([peptide,str(start),str(var),hom_origins[ind],'homozygous']))
    return(variant_info)

def find_str(s, char):
    '''find the position of the peptide inside the protein sequence; return the start position
    '''
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
    parser.add_argument('--report', help='File variant report', required=True)
    parser.add_argument('--pfasta', help='Output protein fasta file', required=False)
    parser.add_argument('--vpep', help='Output variant peptide file', required=False)
    args=vars(parser.parse_args())
    print("Importing data...")
    #import the fasta
    seqs=read_in_fasta(args['fasta'])
    #import the variant information and start codon postition (if applicable)
    var_info=pd.read_csv(args['report'],delim_whitespace=True,header=0,names=['id','haplotype','pos_het','het_origin','pos_hom','hom_origin'])
    print("Preparing data...")
    #merge
    df=pd.merge(var_info,seqs,on=['id','haplotype'])
    #renumber the variant positions to protein position
    for colname in ['pos_het','het_origin','pos_hom','hom_origin']:
        df[colname]=df[colname].str.split(',')
        if colname in ['pos_het','pos_hom']:
            df[colname]=df[colname].apply(renumber_variants)
    print("Predicting protein and peptides...")
    #predict protein sequence
    df['protein_sequence']=df['sequence'].apply(translate)
    #digest protein sequence and get peptides with variants in them
    df['variant_peptides']=df.apply(lambda x: digest(x['protein_sequence'],x['pos_het'],x['het_origin'],x['pos_hom'],x['hom_origin']),axis=1)
    #write information to fasta
    print("Writing to files...")
    if args['pfasta']:
        write_to_fasta(df[['id','haplotype','protein_sequence']],args['pfasta'])
    if args['vpep']:
        var_peps= [item for sublist in list(df['variant_peptides']) for item in sublist]
        f=open(args['vpep'],'w')
        f.writelines(','.join(['peptide','pep_start_in_protein','var_pos_in_protein','variant_origin','var_type'])+'\n')
        f.writelines('\n'.join(var_peps))
