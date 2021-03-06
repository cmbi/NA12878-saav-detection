#!/usr/bin/env python3

'''
this script will recalculate fdr for subsets of peptides, namely variant peptides
variant-free: subset all peptides that were predicted to have variants, recalculate FDR, find the variant peptides back in the new "true positive" set
    - eventually could use cpdt import for this instead
variant-containing: import cpdt files from the variant peptide and decoy script, subset using those, recalculate fdr, see how much I lost

generate qq plots as i go to check if the decoy sets are good or not
'''
import os,sys,re
import pandas as pd
import numpy as np
import file_import
import calculations
import plots
import helper_functions
from collections import Counter
import multiprocessing as mp

def variant_free(df_in,list_var_peps):
    '''
    strategy:
    first recalculate q-values using only peptides predicted by ionbot to have a variant
    then pick out the ones that have genomic evidence that pass new q-value threshold
    '''
    #filter by predicted to have mutation
    df=df_in.loc[df_in["unexpected_modification"].str.contains('[A-Z]->[A-Z]',regex=True)]
    df=calculations.calculate_qvalues(df,decoy_col='DB',score_col='percolator_psm_score')
    indices=np.argwhere(df['q_value']<0.01)
    new_df=df[:int(indices[-1][0])+1] #threshold cut
    out_df=pd.DataFrame()
    for row in new_df.iterrows():
        pep=row[1][3]
        decoy=row[1][6]
        for v in list_var_peps:
            if helper_functions.equivalent_check(pep,v):
                out_df=out_df.append(row[1],ignore_index=True)
    return(out_df,df)

def variant_containing(df_in,cpdt_var,cpdt_rev_var):
    '''
    strategy:
    first find hits that correspond to peptides of interest (variant target and variant decoy variant peptides)
    then re-estimate q-values and take all above threshold

    input: all hits in df, variant peptides and reverse variant peptides
    output: df of all hits that passed threshold
    '''
    # df_in['matched_peptide']=[re.sub('L','I',str(x)) for x in df_in['matched_peptide']]
    # df=df_in.loc[(df_in['matched_peptide'].isin(cpdt_var)) | (df_in['matched_peptide'].isin(cpdt_rev_rp))] #can't do this because I and L need to be equivalent
    # start off the dataframe with the detected variant peptides
    print(len(cpdt_var)) #check the length of the list of true variant peptides
    df=df_in.loc[(df_in['matched_peptide'].isin(cpdt_var))]
    print(df.shape[0]) # should equal the amount of variants found in a particular dataset from variant containing library results
    # df=pd.DataFrame()
    # now add the decoys to the dataframe
    df_decoy=df_in.loc[(df_in['DB']==True)]
    for row in df_decoy.iterrows():
        pep=row[1][3]
        # decoy=row[1][6]
        prot_ids=row[1][9]
        # cpdt=cpdt_rev_var
        cpdt=[]
        if '||' in prot_ids:
            ids=prot_ids.split('||')
        else:
            ids=[prot_ids]
        for i in ids:
            i=helper_functions.get_id(i)
            if i in cpdt_rev_var:
                cpdt+=cpdt_rev_var[i]
        for v in cpdt:
            if pep in v or helper_functions.equivalent_check(pep,v):
                df=df.append(df_in.loc[[row[0]]])
                # break
    df=calculations.calculate_qvalues(df,decoy_col='DB',score_col='percolator_psm_score')
    df.to_csv('test.csv')
    indices=np.argwhere(df['q_value']<0.01)
    return(df[:int(indices[-1][0])+1],df)

def worker_task(csvname):
    dfvf=combi_vf.loc[combi_vf["title"]==csvname.split('.')[0]]
    vf,vff=variant_free(dfvf,list(cpdt_var_vf))
    dfvc=combi_vc.loc[combi_vc["title"]==csvname.split('.')[0]]
    vc,vcf=variant_containing(dfvc,cpdt_var_vc,cpdt_rev)
    return([vf,vc,vff,vcf])

def child_initialize(_combi_vc,_combi_vf,_cpdt_var_vc,_cpdt_var_vf,_cpdt_rev):
     global combi_vc, combi_vf, cpdt_var_vc, cpdt_var_vf, cpdt_rev
     combi_vc = _combi_vc
     combi_vf = _combi_vf
     cpdt_var_vc = _cpdt_var_vc
     cpdt_var_vf = _cpdt_var_vf
     cpdt_rev = _cpdt_rev

def main(combi_vf,combi_vc,cpdt_var_vf,cpdt_var_vc,cpdt_rev_var_file):
    '''
    input:
    - dataframes of all ionbot files (vf and vc)
    - counters of observed variant peptides (vf and vc)
    - file of decoys specifically for the variant peptides (vc)
        - for now only with use with vc FDR re-estimation, but may later add an additional input for the reference counterparts for these decoys in a more specific FDR re-estimation for vf
    
    output: new counters that passed FDR threshold
    '''
    # parser = argparse.ArgumentParser(description='fdr recalculation')
    # parser.add_argument('--cpdt', help='cpdt file variant peptides', required=True)
    # parser.add_argument('--rev', help='cpdt file variant reverse peptides', required=True)
    # parser.add_argument('--vc', help='Directory combi variant-free search results', required=True)
    # parser.add_argument('--vf', help='Directory combi variant-containing search results', required=True)
    # args = vars(parser.parse_args()) 

    # directory_= os.fsencode(csvpath)
    cpdt_rev=file_import.import_cpdt_simple(cpdt_rev_var_file) #import all the variant peptide decoys
    rev_list=[item for sublist in cpdt_rev.values() for item in sublist] #unpack them into a simple list
    ## make L and I equivalent
    # cpdt_var_vc_rp=[w.replace('L','I') for w in cpdt_var_vc] #should not be necessary for the peptides that were taken straight from the results
    # cpdt_rev_rp=[w.replace('L','I') for w in rev_list]
    # cpdt=file_import.import_cpdt(args['cpdt'],False)
    vf=pd.DataFrame()
    vc=pd.DataFrame()
    full_vf=pd.DataFrame()
    full_vc=pd.DataFrame()
    pool=mp.Pool(processes=14,initializer=child_initialize, initargs=(combi_vc,combi_vf,cpdt_var_vc,cpdt_var_vf,cpdt_rev))
    results=[pool.apply_async(worker_task,args=(csvname,)) for csvname in combi_vf["title"].unique()]
    output = [p.get() for p in results]
    for o in output:
        vf=pd.concat([vf,o[0]])
        vc=pd.concat([vc,o[1]])
        full_vf=pd.concat([full_vf,o[2]])
        full_vc=pd.concat([full_vc,o[3]])
    plots.plot_target_decoy(full_vf,'variant_free_ppplot.png') #if this plot is bad, the results are also bad
    plots.plot_target_decoy(full_vc,'variant_containing_ppplot.png')
    #then go through the list and make a new counter
    vfc=Counter(dict(vf['matched_peptide'].value_counts()))
    vcc=Counter(dict(vc['matched_peptide'].value_counts()))
    print(str(sum(vcc.values()))+' variant containing and '+str(sum(vfc.values()))+' variant free variants remaining after FDR correction.')
    #print reports of which proteins they matched to - output both peptide and protein id
    out_df=vf[['proteins','matched_peptide']]
    out_df=out_df.append(vc[['proteins','matched_peptide']])
    out_df=out_df.drop_duplicates()
    out_df.to_csv('variant_pep_report.txt',sep='\t',index=False,header=True)
    return(vfc,vcc)
