#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import phylopandas as ph
import argparse

def check_shift(transcript,protein):
    protein_temp=protein[1:]
    for shift in list(range(1, 10)):
        translated=translate(transcript[shift:])
        if translated[:10]==protein_temp[:10]:
            return(shift)
    return(99)

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

def main():
    parser = argparse.ArgumentParser(description='Get shifts for alternative CDS starts')
    parser.add_argument('--tfasta', help='Transcript fasta: gencode pc_transcripts', required=True)
    parser.add_argument('--pfasta', help='Protein fasta: gencode pc_translations', required=True)
    parser.add_argument('--out', help='Output file for phased exons', required=True)
    args=vars(parser.parse_args())
    proteins=ph.read_fasta(args['pfasta'])
    altstart=proteins[proteins['sequence'].str.startswith('X')]
    transcripts=ph.read_fasta(args['tfasta'])
    transcripts['id']=transcripts['id'].str.split('|').apply(lambda x: x[0])
    altstart=pd.merge(altstart[['id','sequence']],transcripts[['id','sequence']],on='id',suffixes=('_pr','_tr'))
    altstart['shift']=altstart.apply(lambda x: check_shift(x['sequence_tr'],x['sequence_pr']),axis=1)
    altstart[['id','shift']].to_csv(args['out'],index=False)

if __name__ == "__main__":
    main()