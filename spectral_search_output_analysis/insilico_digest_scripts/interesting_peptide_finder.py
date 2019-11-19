#!/usr/bin/env python3

import re, collections, sys
import multiprocessing as mp
import pandas as pd
import argparse

def child_initialize(_ref,_custom):
     global ref, custom
     ref = _ref
     custom = _custom
    
def snvfinder(pid):
    try:
        ref_df=ref[ref['protein']==pid]
        assert (ref_df.shape[0]>0), 'not getting reference df'
        custom_df=custom[custom['protein'].str.contains(pid)]
        saav_list=[]
        for custom_id in custom_df['protein'].unique():
            var_df=custom_df[custom_df['protein']==custom_id]
            assert (var_df.shape[0]>0), 'protein not found in variant-containing pep lib'
            variant_peplist=set(list(var_df['sequence']))
            reference_peplist=set(list(ref_df['sequence']))
            intersect=variant_peplist.intersection(reference_peplist)
            if len(intersect)!=0: # should have some peptides in common
                dif=variant_peplist.difference(reference_peplist) #save all the peptides that are in the custom but not the reference for this particular protein
                for pep in dif:
                    if pep not in ref['sequence']: #need to check if variant peptide is not already present in the reference dictionary
                        start=list(var_df.loc[var_df['sequence']==pep,'start']) #fetch start position (not probability)
                        assert (all_same(start)), start
                        # if float(prob)>=0.05: #higher than cutoff probability
                        isSNV,counterpart,tup=determine_snv(pep,reference_peplist)
                        if isSNV: #counterpart not filtered (optionally can filter that it does not show up in another gene)
                            saav_list.append(','.join([custom_id,pep,counterpart,tup,str(start[0])]))
    except Exception as e:
        raise e
    return('\n'.join(saav_list))

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
    mismatch=[]
    il=False #TO IMPLEMENT: get statistics about how many times we have an IL substitution as the only substitution
    for peptide_comp in plist:
        if len(peptide_comp)==len(peptide_query) and peptide_comp!=peptide_query:
            for a, b in zip(peptide_query,peptide_comp):
                if a!=b and a!='*' and b!='*':
                    if not (a=='I' and b=='L') and not (a=='L' and b=='I'):
                        mismatch.append(f"{a},{b}")
                    else:
                        il=True
            if len(mismatch)==1:
                return (True,peptide_comp,mismatch[0])
    return(False,'','')

def read_csv(csvfile):
    df=pd.read_csv(csvfile)#,names=['protein','sequence','start']
    df['protein']=df['protein'].apply(get_id)
    return(df)

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
    parser.add_argument('--ref', help='Reference csv file', required=True)
    parser.add_argument('--var', help='Variant containing csv file', required=True)
    parser.add_argument('--out', help='Output file snv differing', required=False)
    parser.add_argument('--debug', help='debug option',action='store_true', required=False)
    args=vars(parser.parse_args())
    print('Reading in files...')
    ref=read_csv(args['ref'])
    custom=read_csv(args['var'])
    # if args['debug']:
        # ref=ref.head(100)
    print('Searching for variant peptides...')
    pool=mp.Pool(processes=int(args['cpu']),initializer=child_initialize, initargs=(ref,custom,))
    results=[pool.apply_async(snvfinder,args=(pid,)) for pid in ref['protein'].unique()]
    output = [p.get() for p in results]
    print('Writing to file...')
    if args['out']:
        write_csv(output,args['out'])