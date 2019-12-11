#!/usr/bin/env python3

import re, collections, sys
import multiprocessing as mp
import pandas as pd
import argparse

def child_initialize(_ref,_custom):
     global ref, custom
     ref = _ref
     custom = _custom
    
# def snvfinder(pid):
#     try:
#         ref_df=ref[ref['protein']==pid]
#         assert (ref_df.shape[0]>0), 'not getting reference df'
#         custom_df=custom[custom['protein'].str.contains(pid)]
#         saav_list=[]
#         for custom_id in custom_df['protein'].unique():
#             var_df=custom_df[custom_df['protein']==custom_id]
#             assert (var_df.shape[0]>0), 'protein not found in variant-containing pep lib'
#             variant_peplist=set(list(var_df['sequence']))
#             reference_peplist=set(list(ref_df['sequence']))
#             intersect=variant_peplist.intersection(reference_peplist)
#             if len(intersect)!=0: # should have some peptides in common
#                 dif=variant_peplist.difference(reference_peplist) #save all the peptides that are in the custom but not the reference for this particular protein
#                 for pep in dif:
#                     if pep not in ref['sequence']: #need to check if variant peptide is not already present in the reference dictionary
#                         start=list(var_df.loc[var_df['sequence']==pep,'start']) #fetch start position (not probability)
#                         assert (all_same(start)), start
#                         # if float(prob)>=0.05: #higher than cutoff probability
#                         isSNV,counterpart,tup=determine_snv(pep,reference_peplist)
#                         if isSNV: #counterpart not filtered (optionally can filter that it does not show up in another gene)
#                             saav_list.append(','.join([custom_id,pep,counterpart,tup,str(start[0])]))
#     except Exception as e:
#         raise e
#     return('\n'.join(saav_list))

def snvfinder(var_peptides,ref_peptides):
    try:
        assert len(var_peptides.intersection(ref_peptides))>0, var_sequence #make sure there is some overlap between the sequences
        final_variant_peptides=[]
        for pep in var_peptides:
            counterpart=determine_snv(pep,ref_peptides)
            if counterpart!='':
                final_variant_peptides.append(result)
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
    if df['description'].str.contains('hap:').all():
        df['haplotype']=df['description'].str.split(' hap:').apply(lambda x: x[1])
        return(df[['id','haplotype','sequence']])
    return(df[['id','sequence']])

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
    parser.add_argument('--cpu',help='CPU number',required=True)
    parser.add_argument('--ref', help='Reference protein fasta', required=True)
    parser.add_argument('--var', help='Variant containing protein fasta', required=True)
    parser.add_argument('--prelim', help='Preliminary variant peptides generated by "translate_annotated"', required=False)
    parser.add_argument('--out', help='Output file snv differing', required=False)
    parser.add_argument('--debug', help='debug option',action='store_true', required=False)
    args=vars(parser.parse_args())
    print('Reading in files...')
    ref=read_fasta(args['ref'])
    custom=read_fasta(args['var'])
    merged=pd.merge(custom,ref,on='id',suffixes=('_vc','_vf'))
    print('Searching for variant peptides...')
    merged['variant_peps']=merged.apply(lambda x: snvfinder(x['peptides_vc'],x['peptides_vf']),axis=1)
    df=merged[['id','haplotype','variant_peps']]
    df=df[df['variant_peps'].map(lambda d: len(d)) > 0] #remove proteins without any variant peptides
    #seperate variant peptides into their own lines
    unstacked=df.apply(lambda x: pd.Series(x['variant_peps']),axis=1).stack().str.strip().reset_index(level=1, drop=True)
    unstacked.name='peptide'
    df=df.drop('variant_peps',axis=1).join(unstacked)
    #separate the peptide and the substitution
    df['peptide']=df['peptide'].str.split('|')
    df[['peptide','substitution']]=pd.DataFrame(df['peptide'].values.tolist(),index=df.index)
    if args['prelim']:
        origins=pd.read_csv(args['prelim'])
        df=pd.merge(df,origins,on='peptide')
    # if args['debug']:
        # ref=ref.head(100)
    # pool=mp.Pool(processes=int(args['cpu']),initializer=child_initialize, initargs=(ref,custom,))
    # results=[pool.apply_async(snvfinder,args=(pid,)) for pid in ref['protein'].unique()]
    # output = [p.get() for p in results]
    print('Writing to file...')
    if args['out']:
        write_csv(df,args['out'])