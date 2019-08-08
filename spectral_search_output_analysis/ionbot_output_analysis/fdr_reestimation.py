#!/usr/bin/env python3

'''
this script will recalculate fdr for subsets of peptides, namely variant peptides
variant-free: subset all peptides that were predicted to have variants, recalculate FDR, find the variant peptides back in the new "true positive" set
    - eventually could use cpdt import for this instead
variant-containing: import cpdt files from the variant peptide and decoy script, subset using those, recalculate fdr, see how much I lost

generate qq plots as i go to check if the decoy sets are good or not
'''
import pandas as pd
import numpy as np
import file_import
import calculations
import plots

def variant_free(df_in,list_var_peps):
    '''
    strategy:
    first recalculate q-values using only peptides predicted by ionbot to have a variant
    then pick out the ones that have genomic evidence that pass new q-value threshold
    '''
    #filter by predicted to have mutation
    df=df_in.loc[df_in["unexpected_modification"].str.contains('[A-Z]->[A-Z]',regex=True)]
    df=calculations.calculate_qvalues(df,decoy_col='DB',score_col='percolator_psm_score')
    plots.plot_target_decoy(df,'variant_free_ppplot.png') #if this plot is bad, the results are also bad
    indices=np.argwhere(df['q_value']<0.01)
    new_df=df[:int(indices[-1][0])+1] #threshold cut
    out_df=pd.DataFrame()
    for row in new_df.iterrows():
        pep=row[1][3]
        decoy=row[1][6]
        for v in list_var_peps:
            if helper_functions.equivalence_check(pep,v):
                out_df=out_df.append(row[1],ignore_index=True)
    return(out_df)

def variant_containing(df_in,cpdt_var,cpdt_rev_var):
    '''
    strategy:
    first find hits that correspond to peptides of interest (variant and decoy variant peptides)
    then re-estimate q-values and take all above threshold

    input: all hits in df, variant peptides and reverse variant peptides
    output: df of all hits that passed threshold
    '''
    # df=df_in.loc[(df_in['matched_peptide'].isin(cpdt_var.values())) | (df_in['matched_peptide'].isin(cpdt_rev_var.values()))] #can't do this because I and L need to be equivalent
    df=pd.DataFrame()
    for row in df_in.iterrows():
        pep=row[1][3]
        decoy=row[1][6]
        if decoy=='D':
            cpdt=cpdt_rev_var
        else:
            cpdt=cpdt_var
        for v in cpdt.values():
            if helper_functions.equivalence_check(pep,v):
                df=df.append(row[1],ignore_index=True)
    df=calculations.calculate_qvalues(df,decoy_col='DB',score_col='percolator_psm_score')
    plots.plot_target_decoy(df,'variant_containing_ppplot.png')
    indices=np.argwhere(df['q_value']<0.01)
    return(df[:int(indices[-1][0])+1])


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
    cpdt_rev=file_import.import_cpdt(cpdt_rev_var_file,False)
    # cpdt=file_import.import_cpdt(args['cpdt'],False)
    vf=pd.DataFrame()
    vc=pd.DataFrame()
    for csvname in combi_vf["title"].unique(): #calc FDR in a dataset-specific manner
        df=combi_vf.loc[combi_vf["title"]==csvname.split('.')[0]]
        vf=pd.concat([vf,variant_free(df,list(cpdt_var_vf))])
        vc=pd.concat([vc,variant_containing(df,cpdt_var_vc,cpdt_rev)])
    #then go through the list and make a new counter
    vfc=Counter(dict(vf['matched_peptide'].value_counts()))
    vcc=Counter(dict(vc['matched_peptide'].value_counts()))
    return(vfc,vcc)
