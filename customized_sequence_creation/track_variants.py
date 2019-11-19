#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import gffpandas.gffpandas as gffpd
import argparse

'''
this script will use as input:
    - the report output from combine_exons
    - CDS
    - file with identifiers of the proteins and peptides from fdr_reestimation

will return the SNP information (first file) from those variants that were found in the peptides (last file)

1. read in CDS from ANGEL output (suffix "final.cds")
2. match CDS with transcript sequence on record
3. renumber the variant position to CDS
4. renumber the variant position to protein
5. result: peptide - SNP linked

'''

def read_filter_cds(cdsfile):
    if cdsfile[-3:]=='gff3':
        df=gffpd.read_csv(cdsfile)
    else:
        df=pd.read_csv(cdsfile,sep=' ')
    return(cds)



def map_seq(reportdf):
    '''
    separate the list of proteins, so that it is 1 protein, 1 peptide format
    '''

def gather_prots(list_prots):
    out=set()
    for p in list_prots:
        if '||' in p:
            s=set(p.split('||'))
            out=out.union(s)
        else:
            out.add(p)
    return out

def protein_filter(prots):
    o=set()
    for p in prots:
        o.add(get_id(p))
    return(o)

def get_id(idstring):
    i=idstring
    if '|m.' in i:
        i=i.split('|m.')[0]
    elif 'ENSP' in i:
        i=i.split('|')[1]
    return(i)

def main():
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--var', help='', required=True)
    parser.add_argument('--cds', help='angel final cds file for all the CDS that are not in ref', required=True)
    args=parser.parse_args()

    cds_ont=read_filter_cds(args.cdsont,prots) #protein \t cds
    cds_ref=read_filter_cds(args.cdsref,prots)
    report=pd.read_csv(args.var,sep='\t',columns=["transcript","local_pos","chrom_pos"]) #transcript \t var positions on transcript \t chr|chrom_positions
    transcripts=protein_filter(prots)
    report=report.loc[report["transcript"].isin(transcripts)]
    
    