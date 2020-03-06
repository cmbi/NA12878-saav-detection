#!/usr/bin/env python3

import re, collections, sys
import multiprocessing as mp
import pandas as pd
import argparse
import phylopandas as ph
from pyteomics.parser import cleave

def snvfinder(var_peptides,ref_peptides,tid):
    try:
        final_variant_peptides=[]
        for pep in var_peptides:
            counterpart=determine_snv(pep,ref_peptides)
            if counterpart!='':
                final_variant_peptides.append(f"{pep}|{counterpart}")
    except Exception as e:
        raise e
    return(final_variant_peptides)

def all_same(items):
    return(all(x == items[0] for x in items))

def get_id(idstring):
    if '|m.' in idstring:
        outstring=idstring.split('|m.')[0]
    elif 'ENSP' in idstring:
        tid=idstring.split('|')[1]
        if 'Random' in idstring:
            prefix=idstring.split('_')[0]
            outstring=prefix+'_'+tid
        else:
            outstring=tid
    else:
        outstring=idstring.strip()
    return(outstring)

def determine_snv(peptide_query,plist):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: boolean, the reference peptide, and the 
    '''
    il=False #TO IMPLEMENT: get statistics about how many times we have an IL substitution as the only substitution
    for peptide_comp in plist:
        if len(peptide_comp)==len(peptide_query) and peptide_comp!=peptide_query:
            mismatch=[]
            for a, b in zip(peptide_query,peptide_comp):
                if a!=b and a!='*' and b!='*':
                    if not (a=='I' and b=='L') and not (a=='L' and b=='I'):
                        mismatch.append(f"{a},{b}")
                    else:
                        il=True
            if len(mismatch)==1:
                return (f'{peptide_comp}|{mismatch[0]}')
    return('')

def digest(protein_sequence):
    '''digests the protein sequences with trypsin
    '''
    seq_cut = cleave(protein_sequence, '[KR]', 2)
    plist=set() #to prevent duplicates
    for peptide in seq_cut:
        if peptide in plist or len(peptide) < 6:
            continue
        plist.add(peptide)
    return(plist)

def read_fasta(fafile):
    df=ph.read_fasta(fafile)#,names=['protein','sequence','start']
    df['id']=df['id'].apply(get_id)
    df['peptides']=df['sequence'].apply(digest)
    if df['id'].str.contains('\|h.').all():
        df['haplotype']=df['id'].str.split('\|h.')#.apply(lambda x: x[-1])
        df[['id','haplotype']]=pd.DataFrame(df.haplotype.values.tolist(), index= df.index)
        return(df[['id','haplotype','sequence','peptides']])
    return(df[['id','sequence','peptides']])

def write_csv(saav_list,outfile):
    '''
    input format: {id:{peptide:probability}}
    output: new csv
    '''
    f=open(outfile,'w')
    f.writelines(','.join(['protein', 'variant','counterpart', 'start'])+'\n')
    for s in saav_list:
        if s!='':
            f.writelines(s+'\n')
    return(f'csv file written to {outfile}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='track variants')
    # parser.add_argument('--cpu',help='CPU number',required=True)
    parser.add_argument('--ref', help='Reference protein fasta', required=True)
    parser.add_argument('--var', help='Variant containing protein fasta', required=True)
    parser.add_argument('--prelim_gen', help='Preliminary variant peptides generated by "translate_annotated"', required=False)
    parser.add_argument('--prelim_ont', help='Preliminary variant peptides generated by "translate_annotated"', required=False)
    parser.add_argument('--out', help='Output file snv differing', required=False)
    parser.add_argument('--debug', help='debug option',action='store_true', required=False)
    args=vars(parser.parse_args())
    print('Reading in files...')
    ref=read_fasta(args['ref'])
    custom=read_fasta(args['var'])
    merged=pd.merge(custom,ref,on='id',suffixes=('_vc','_vf'))
    if args['debug']:
        merged=merged.head(500)
    print('Searching for variant peptides...')
    merged['variant_peps']=merged.apply(lambda x: snvfinder(x['peptides_vc'],x['peptides_vf'],x['id']),axis=1)
    df=merged[['id','haplotype','variant_peps']]
    df=df[df['variant_peps'].map(lambda d: len(d)) > 0] #remove proteins without any variant peptides
    #seperate variant peptides into their own lines
    unstacked=df.apply(lambda x: pd.Series(x['variant_peps']),axis=1).stack().str.strip().reset_index(level=1, drop=True)
    unstacked.name='peptide'
    df=df.drop('variant_peps',axis=1).join(unstacked).reset_index(drop=True)
    #separate the peptide and the substitution
    df['peptide']=df['peptide'].str.split('|')
    df[['peptide','ref_counterpart','substitution']]=pd.DataFrame(df['peptide'].values.tolist(),index=df.index)
    #add origin information if available
    if args['prelim_gen']:
        origins=pd.concat([pd.read_csv(args['prelim_gen']), pd.read_csv(args['prelim_ont'])], ignore_index=True)
        df=pd.merge(df,origins[['peptide','variant_origin','var_type']].drop_duplicates(),on='peptide')
    print('Writing to file...')
    if args['out']:
        df.to_csv(args['out'],index=False)
        # write_csv(df,args['out'])