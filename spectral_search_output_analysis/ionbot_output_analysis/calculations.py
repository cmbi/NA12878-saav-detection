#!/usr/bin/env python3

import re, os,sys, itertools
import numpy as np
from scipy import stats
import pandas as pd
import helper_functions
import plots

def calculate_qvalues(df_in, decoy_col='decoy', score_col='score',
                      lower_score_is_better=False):
    """
    Calculate q-values using the Target-Decoy Approach.
​
    Positional arguments:
    df_in -- pandas.DataFrame object containing at least the following columns:
    1. Search engine scores, as int or float
    2. Booleans for every search engine score: True if the score corresponds to
    a decoy, False if the score corresponds to a target.
​
    Keyword arguments:
    score_col -- String with the column name for the search engine scores
    (default: 'score')
    decoy_col -- String with the column name for the decoy booleans
    (default: 'decoy')
    lower_score_is_better -- True if a lower search engine score is better.
    (default: False)
​
    More info: Aggarwal,S. and Yadav,A. (2015) In Statistical Analysis in
    Proteomics. https://doi.org/10.1007/978-1-4939-3106-4_7, Chapter 3.1,
    Formulae 1 and 3.
​
    Returns:
    A Pandas Series with the calculated q-values, sorted by the original index.
    """
    df = df_in.sort_values(score_col, ascending=lower_score_is_better)
    df.reset_index(drop=True)
    # if not df[decoy_col].isin(['D']).empty:
    #     df=df.reset_index(drop=True)
    #     df=df.drop(df[(df[decoy_col]!='D') & (df[decoy_col]!='T')].index)
    #     df[decoy_col]=df[decoy_col].map({'D':True,'T':False})
    df['decoy_cumsum'] = df[decoy_col].cumsum()
    df['target_cumsum'] = (~df[decoy_col]).cumsum()
    df['q_value'] = (df['decoy_cumsum'] / df['target_cumsum'])
    return df

def fdr_recalc_variantpep(target,decoy,filename):
    '''filter df by new q-value
    '''
    # print(str(target.shape[0])+' target variant peptides and '+str(decoy.shape[0])+' decoy variant peptides found.')
    df_in=pd.concat([target,decoy],ignore_index=True)
    plots.plot_target_decoy(df_in,filename) #if this plot is bad, the results are also bad
    df=calculate_qvalues(df_in,decoy_col='DB',score_col='percolator_psm_score')
    indices=np.argwhere(df['q_value']<0.01)
    return(df[:int(indices[-1][0])+1]) #threshold cut


def coverage_measure(cpdt_pep,full_seqs):
    #high_cov_vert={}
    high_cov_hor={}
    perc_cov_dist=[]
    vert_cov=[]
    for p,peps in cpdt_pep.items():
        seq=full_seqs[p]
        remains=seq
        count_pep=0
        if len(peps.keys())>0:
            for s,c in peps.items():
                if c>5: #what constitutes a "true" hit
                    count_pep+=c
                    if s in remains:
                        remains=remains.replace(s,'')
                    else:
                        prefix=re.split('R|K',s)
                        for p in prefix:
                            if len(p)>3 and p in remains:
                                remains=remains.replace(p,'')
            perc_cov=float((len(seq)-len(remains))/len(seq))*100
            perc_cov_dist.append(perc_cov)
            vert=float(count_pep/len(peps.keys()))
            if vert>0:
                vert_cov.append(vert)
            if perc_cov>50:
                high_cov_hor[p]=peps
            # if count_pep>100:
            #     high_cov_vert[p]=peps
    return(high_cov_hor,vert_cov,perc_cov_dist)

def calc_nsaf_standard(cpdt_pep,fullseqs):
    '''after cpdt_pep has been filled, find the sum nsaf in order to standardize the abundance scores in calc_mut_abundances'''
    nsaf=0
    for prot,peps in cpdt_pep.items():
        count_peps=0
        for p,ct in peps.items():
            count_peps+=ct
        length_prot=len(fullseqs[prot])
        nsaf+=float(count_peps/length_prot)
    return(nsaf)

def calc_nsaf_protein(pepdict_singleprot,lenprot,sumnsaf):
    sum_nonmut=0
    for normpep,normct in pepdict_singleprot.items():
        sum_nonmut+=normct
    nsaf=float(float(sum_nonmut/lenprot)/sumnsaf)
    return(nsaf)

def calc_mut_abundances(mutant_cpdtpep,cpdtpep,fullseqs):
    '''this function calculates the abundance of the protein and returns that along with the count of the (non)mutant peptide'''
    mut_pep_abundance=[]
    nonmut_pep_abundance=[]
    #nr_mutant=[]
    mut_proteins_detected=set()
    num_peptides=0
    num_occurences=0
    sumnsaf=calc_nsaf_standard(cpdtpep,fullseqs)
    for prot,peps in mutant_cpdtpep.items():
        if '_h' in prot:
            stem=prot.split('_h')[0]
        else:
            stem=prot
        if stem in cpdtpep:
            nsafnonmut=calc_nsaf_protein(cpdtpep[stem],len(fullseqs[stem]),sumnsaf)
            sum_mut=0 #total number of detected mutant peptides
            sum_nonmut=0
            for pep,ct in peps.items():
                sum_mut+=ct
                if ct>0:
                    mut_proteins_detected.add(prot)
                    num_occurences+=ct
                    num_peptides+=1
                    # mut_pep_abundance.append((nsafnonmut,ct)) #uncomment to calculate frequencies per peptide instead of per protein
            #calculate count of non-mutant peptides
            for normpep,normct in cpdtpep[stem].items():
                sum_nonmut+=normct
            lennonmut=len(fullseqs[stem])
            for normpep,normct in cpdtpep[stem].items():
                sum_nonmut+=normct
                # nonmut_pep_abundance.append((nsafnonmut,normct)) #uncomment to calculate frequencies per peptide instead of per protein
            # nsafnonmut=float(float(sum_nonmut/lennonmut)/sumnsaf)
            nonmut_pep_abundance.append((nsafnonmut,sum_nonmut))
            if sum_mut>0: #only record the proteins with at least 1 detected variant peptide
                #nr_mutant.append(sum_mut)
                mut_pep_abundance.append((nsafnonmut,sum_mut))
                #mut_pep_abundance.append(sum_mut+sum_nonmut)
    print("Total of "+str(num_occurences)+" occurances of "+str(num_peptides)+" peptides from "+str(len(mut_proteins_detected))+" proteins were detected")
    return(mut_pep_abundance,nonmut_pep_abundance)

def calc_pep_counts(mutant_cpdtpep,counterpart_cpdtpep,mutant_probs):
    '''directly measure the observed counts of the saav peptides and their reference counterparts
    for later implementation: coloring similar to the direct comparison discrepant peptide graph, with colors corresponding to the abundance of peptides
    '''
    counts=[]
    probs=[]
    observed_subs=helper_functions.initiate_counter()
    for prot,peps in mutant_cpdtpep.items():
        if prot in counterpart_cpdtpep:
            peps_cpt=counterpart_cpdtpep[prot]
            for pep in peps:
                cpt_pep,sub=helper_functions.determine_snv(pep,peps_cpt)
                if cpt_pep in peps_cpt:
                    tuptoadd=(peps[pep],peps_cpt[cpt_pep])
                    probtup=(mutant_probs[pep],peps[pep])
                    if tuptoadd[0]!=0: #only if at least 1 variant peptide detected
                        counts.append(tuptoadd)
                        probs.append(probtup)
                        observed_subs[sub]+=1
    # df,dfall=helper_functions.counter_to_df(observed_subs)
    return(counts,probs,observed_subs)

def saav_counts(variantdf,counterpartdf,peptide_colname='peptide',observed=False):
    '''input the variant and counterpart dataframes to get a count of which/how many AA substitutions occur'''
    all_subs=helper_functions.initiate_counter()
    count_subs=[]
    for protname in variantdf['proteins'].unique():
        slice_var=variantdf[variantdf['proteins']==protname]
        slice_ctp=counterpartdf[counterpartdf['proteins']==protname]
        variant=slice_var[peptide_colname].unique()
        counterpart=slice_ctp[peptide_colname].unique()
        for var in variant:
            cpt_pep,sub=helper_functions.determine_snv(var,counterpart)
            if sub!='':
                all_subs[sub]+=1
                if observed:
                    count_var=slice_var[slice_var[peptide_colname]==var].shape[0] #how many of this variant peptide was detected
                    count_ctp=slice_ctp[slice_ctp[peptide_colname]==cpt_pep].shape[0] #how many of the counterpart was detected
                    count_subs.append((count_var,count_ctp))
    if observed:
        return(all_subs,count_subs)
    # for prot,peps in mutant_cpdtpep.items():
    #     if prot in counterpart_cpdtpep:
    #         peps_cpt=counterpart_cpdtpep[prot]
    #         for pep in peps:
                
    # serall=pd.Series(list(all_subs.values()),index=pd.MultiIndex.from_tuples(all_subs.keys()))
    # dfall=serall.unstack()
    return(all_subs)

def calculate_correlation(mut_pep_abundance,cpt_abundance,nonmut_pep_abundance):
    dt=np.dtype('float,int')
    variant = np.array(mut_pep_abundance,dtype=dt)
    nonvariant=np.array(nonmut_pep_abundance,dtype=dt)
    cor_var= np.corrcoef(variant['f0'],variant['f1'])
    cor_cpt=np.corrcoef(variant['f0'],variant['f1'])
    cor_nonvar= np.corrcoef(nonvariant['f0'],nonvariant['f1'])
    # cor_var_nonvar=np.corrcoef(variant['f0'],nonvariant['f0'])
    print("correlation between variant peptides abundance and total protein abundance: "+str(cor_var[1][0]))
    print("correlation between counterpart peptide abundance and total protein abundance:"+str(cor_cpt[1][0]))
    print("correlation between non-variant peptide abundance and total protein abundance: "+str(cor_nonvar[1][0]))
    return(0)


def r2(x, y):
    # return(stats.pearsonr(x, y)[0] ** 2)
    tup=stats.pearsonr(x, y)
    rsq=tup[0] ** 2
    return(rsq,tup[1])