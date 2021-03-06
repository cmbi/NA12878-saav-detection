#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import gffpandas.gffpandas as gffpd
import argparse

def get_ends(tstarts,blocksizes):
    ends=[]
    for i,t in enumerate(tstarts):
        ends.append(int(t)+int(blocksizes[i]))
    return(ends)

def get_position(cds_start,blocksizes):
    '''return the cds that has the start codon 
    and return the shift between the start codon and beginning of corresponding cds'''
    current_blocksize=0
    for index,blocksize in enumerate(blocksizes):
        if cds_start > current_blocksize and cds_start <= current_blocksize + int(blocksize):
            shift_from_end_exon=current_blocksize+int(blocksize)-cds_start+1
            return(index,(int(blocksize)-shift_from_end_exon))
        current_blocksize+=int(blocksize)

def update_blocks(cds_start,blocksizes,genome_starts,genome_ends,strand):
    ''' should return 3 lists
    list 1: updated block sizes
    list 2: updated start positions
    list 3: consecutive numbers in a list the same length as the other 2 lists (for ranking)
    '''
    exon_num,shift=get_position(cds_start,blocksizes)
    genome_starts=genome_starts[exon_num:]
    genome_ends=genome_ends[exon_num:]
    if strand=='+':
        genome_starts[0]=int(genome_starts[0])+shift
    elif strand=='-':
        genome_ends[0]=int(genome_ends[0])-shift
    exonnumbers=list(range(1, len(genome_starts)+1))
    return(genome_starts,genome_ends,exonnumbers)

def convert_to_gff(df,gffoutfile):
    new_df=pd.DataFrame()
    random_df=pd.DataFrame()
    gff=gffpd.Gff3DataFrame()
    new_df['seq_id']=df[13]
    new_df['transcript']=df[9]
    new_df['score']='.'
    new_df['strand']=df[8]
    new_df['phase']='.'
    random_df['transcript']=df[9]
    random_df['start']=df['tStarts'].apply(lambda x: [int(d)+1 for d in x]) #0- to 1-based coordinates
    random_df['end']=df['tEnds']
    random_df['exon_number']=df['exon_number']
    new_df['source']='ONT'
    new_df['type']='CDS'
    start_ends=random_df.set_index(['transcript']).apply(pd.Series.explode).reset_index() #separate so one start/end per line
    new_df=new_df.merge(start_ends, on='transcript').reset_index(drop=True)
    new_df['attributes']='transcript_id='+new_df['transcript']+';exon_number='+new_df['exon_number'].astype(str)+';'
    new_df['end']=new_df['end'].astype(int)
    gff.df=new_df[['seq_id','source','type','start','end','score','strand','phase','attributes']]
    gff.header=f'##gff-version 3\n#description: CDS regions of sample genome\n'
    gff.to_gff3(gffoutfile)

def main():
    parser = argparse.ArgumentParser(description='PSL to (CDS) gff3 converter')
    parser.add_argument('--sq', help='SQANTI2 output classification file (for CDS information)', required=True)
    parser.add_argument('--psl', help='FLAIR output PSL file', required=True)
    parser.add_argument('--out', help='Output gff3', required=True)
    args=vars(parser.parse_args())
    classification=pd.read_csv(args['sq'],sep='\t')
    # classification=classification[(classification['coding']=='coding')&(~classification['subcategory'].str.contains('fragment'))] #shouldn't need this; filtered in the input
    classification=classification[['isoform','CDS_start']]
    psl=pd.read_csv(args['psl'],sep='\t',header=None)
    psl['isoform']=psl[9]
    psl=psl.merge(classification,on='isoform') #check here what I am getting rid of, should not be less than classification.shape[0]
    psl['tStarts']=psl[20].str.split(',').apply(lambda x: list(filter(None, x)))
    psl['tStarts']=psl.apply(lambda x: x['tStarts'][::-1] if x[8]=='-' else x['tStarts'],axis=1)
    psl['blocksize']=psl[18].str.split(',').apply(lambda x: list(filter(None, x)))
    psl['blocksize']=psl.apply(lambda x: x['blocksize'][::-1] if x[8]=='-' else x['blocksize'],axis=1)
    psl['tEnds']=psl.apply(lambda x: get_ends(x['tStarts'],x['blocksize']),axis=1)
    psl['exon_number']=psl.apply(lambda x: update_blocks(x['CDS_start'],x['blocksize'],x['tStarts'],x['tEnds'],x[8]),axis=1)
    psl[['tStarts','tEnds','exon_number']]=pd.DataFrame(psl['exon_number'].tolist(),index=psl.index)
    convert_to_gff(psl,args['out'])

if __name__ == "__main__":
    main()