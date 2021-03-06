#!/usr/bin/env python3

import re, os,sys, itertools, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
import plots
import helper_functions
import calculations
import file_import
import logging
# log = logging.getLogger(__name__)

def protein_support(df_vf,df_vc):
    '''calculate support for proteins
    input the "proteins" column of ionbot only
    plot the protein support
    '''
    support_vf=Counter([x for l in df_vf.tolist() for x in l.split('||')])
    unamb_support_vf=Counter(df_vf[~df_vf.str.contains('\|{2}',regex=True)].tolist())

    support_vc=Counter([x for l in df_vc.tolist() for x in l.split('||')])
    unamb_support_vc=Counter(df_vc[~df_vc.str.contains('\|{2}',regex=True)].tolist())

    plots.plot_support(support_vf,unamb_support_vf,'protein_evidence_varfree.png')
    plots.plot_support(support_vc,unamb_support_vc,'protein_evidence_varcont.png')

def origin_info_fetch(df_vf,df_vc,gff,bed):
    '''
    '''
    def fetch_umambig(prot_string):
        if '||' not in prot_string:
            if '|h' in prot_string:
                return(prot_string.split('|h')[0])
            else:
                return(prot_string.split('((')[0])
        return('ambiguous')
    df_vf['transcript_id']=df_vf['proteins'].apply(fetch_umambig)#lambda x: x.split('|h')[0] if '|h' in x else x.split('((')[0])
    df_vc['transcript_id']=df_vc['proteins'].apply(fetch_umambig)#lambda x: x.split('|h')[0] if '|h' in x else x.split('((')[0])
    origin_info=file_import.create_chromosome_reference(gff,bed)
    vf=df_vf.merge(origin_info,on='transcript_id').drop_duplicates(subset=('title'))
    vc=df_vc.merge(origin_info,on='transcript_id').drop_duplicates(subset=('title'))
    strand_vf=Counter(dict(vf['strand'].value_counts(normalize=True)))
    strand_vc=Counter(dict(vc['strand'].value_counts(normalize=True)))
    chrom_vf=Counter(dict(vf['chromosome'].value_counts(normalize=True)))
    chrom_vc=Counter(dict(vc['chromosome'].value_counts(normalize=True)))
    vf[vf['chromosome'].isin(['chr12','chr14','chr22'])].to_csv('different_chr_origins_vf.csv')
    vc[vc['chromosome'].isin(['chr12','chr14','chr22'])].to_csv('different_chr_origins_vc.csv')
    plots.plot_chromosomal_dist(chrom_vc,chrom_vf)
    plots.plot_strand_dist(strand_vc,strand_vf)

def dict_source_bin(df_vf,df_vc,figname_vf,figname_vc):
    vf=Counter(dict(df_vf.value_counts()))
    vc=Counter(dict(df_vc.value_counts()))
    df_vf.value_counts().to_csv(figname_vf+'.csv')
    df_vc.value_counts().to_csv(figname_vc+'.csv')
    plots.plot_source_piechart(vf,figname_vf)
    plots.plot_source_piechart(vc,figname_vc)        

def discrepancy_check(vf,vc,nonvar_vf,nonvar_vc,var_peps):
    '''check out the differences in identifications between the 2 combination dictionaries
    why doesn't ionbot catch everything? look at the ones that it does not catch but the variant-containing dictionary does
    plot lengths of the missed/caught peptides (longer than average?)
    plot unexpected modifications of the missed peptides (more unexpected modifications than average?)'''
    #abbreviate the results, group by matched peptide
    # aggregation_rules={'percolator_psm_score_best':'mean','unexpected_modification':helper_functions.longest,'proteins':helper_functions.longest,'matched_peptide':helper_functions.longest}
    # ibdf_vf=vf.groupby('peptide').aggregate(aggregation_rules)
    # ibdf_vc=vc.groupby('peptide').aggregate(aggregation_rules)
    #convert everything to sets and find the variants that are/are not in common
    merged=pd.merge(vf,vc, on='title', how='outer', suffixes=('_vf','_vc'), indicator=True)
    vf_only=merged[merged['_merge'] == 'left_only']
    vc_only=merged[merged['_merge'] == 'right_only']
    vc_only[['title']].drop_duplicates().merge(vc,on='title').drop_duplicates().merge(nonvar_vf,on='title',suffixes=('_vf','_vc')).drop_duplicates().to_csv('mislabeled_new.csv',index=False)
    overlap=merged[merged['_merge'] == 'both']
    #plot lengths of the variant peptides caught by VC but not VF against nonvariant VC
    vcopl=[len(i) for i in vc_only['matched_peptide_vc'].unique()]
    vcnvpl=[len(i) for i in nonvar_vc['matched_peptide'].unique()]
    vpl=[len(i) for i in var_peps['variant_peptide'].unique()]
    plots.plot_peplengths(helper_functions.normalize_counter(Counter(vcopl)),helper_functions.normalize_counter(Counter(vpl)),helper_functions.normalize_counter(Counter(vcnvpl)))
    # plots.plot_peplengths(helper_functions.normalize_counter(Counter(vc_only['matched_peptide_vc'].str.len().to_list())),helper_functions.normalize_counter(Counter(nonvar_vc['matched_peptide'].str.len().to_list()))) #does not take into account repeating peptides
    #look at unexpected modifications in the mis-labeled VF (since there are no VF exclusive variants)
    vc_unique=pd.merge(vf['title'],vc[['title','variant_peptide','matched_peptide','percolator_psm_score_best','ionbot_score','corr','rt','predicted_tr']], on='title', how='right', indicator=True) #isolate variant containing only
    mislabeled=pd.merge(nonvar_vf[['title','matched_peptide','modifications','unexpected_modification','best_psm','ionbot_score','corr','percolator_psm_score_best','rt','predicted_tr']],vc_unique.loc[vc_unique['_merge']=='right_only',('title','matched_peptide','best_psm','percolator_psm_score_best','ionbot_score','corr','rt_vc')], on='title',suffixes=('_vf','_vc'))
    #check the distributions
    mislabeled.to_csv('mislabeled_varpeps.csv',index=False) #comment this out if not want 
    #mislabeled_small=pd.merge(nonvar_vf[['title','peptide','unexpected_modification']].fillna('none'),vc_unique.loc[vc_unique['_merge']=='right_only','title'], on='title').groupby(['peptide','unexpected_modification']).aggregate(lambda x: x.iloc[0]).fillna("none").reset_index()
    plots.plot_unexpected_mods(mislabeled["unexpected_modification"],nonvar_vf.fillna('').apply(lambda x: helper_functions.categorize_mods(x["unexpected_modification"],x["pred_aa_sub"]),axis=1))
    #direct comparison of scores of peptides
    # rt_obs_df=pd.read_csv(rt_obs)
    #to write to file
    # mislabeled=mislabeled.groupby(['matched_peptide_vc','matched_peptide_vf'])#.aggregate('mean').reset_index() #group by peptide identification
    # mislabeled_distr=mislabeled[["percolator_psm_score_best_vc","percolator_psm_score_best_vf"]].describe()
    #scores_all_filtered=mislabeled.loc[mislabeled["percolator_psm_score_best_vc"]>mislabeled["percolator_psm_score_best_vf"]] #information for scan ids that are higher in variant containing than variant free
    #scores_all_filtered.to_csv('vc_higher_than_vf.csv',index=False) #first print to csv, look at where the vc scores were higher than vf
    plots.plot_ib_scores_directcomp(mislabeled.rename(columns={"matched_peptide_vc":"matched_peptide"})) #direct comparison plot: what scores they had in each of the libraries
    #general comparison of the scores
    plots.plot_ib_scores(vf_only["percolator_psm_score_best_vf"].tolist(),vc_only["percolator_psm_score_best_vc"].tolist(),overlap["percolator_psm_score_best_vc"].tolist(),overlap["percolator_psm_score_best_vf"].tolist(),nonvar_vc["percolator_psm_score_best"].tolist(),nonvar_vf["percolator_psm_score_best"].tolist())
    #to explore: return scan ids and check the ids that were not identified with variant free method in a later function
   
def explore_rt_dist(variant_df,nonvariant_df,plotname):
    variant_df['delta_rt']=abs(variant_df['tr']-variant_df['predicted_tr'])
    nonvariant_df['delta_rt']=abs(nonvariant_df['tr']-nonvariant_df['predicted_tr'])
    #for those variant findings that are in agreement between the 2 methods
    overlap=merged[merged['_merge'] == 'both']
    overlap_agreement=overlap[overlap['peptide_vc']==overlap['peptide_vf']]
    overlap_stats=overlap_agreement['delta_rt'].groupby(['peptide_vc']).describe()
    #for those variant findings that are not in agreement between the 2 methods
    variant_disagreement=overlap[overlap['peptide_vc']!=overlap['peptide_vf']]
    variant_disagreement_stats=variant_disagreement['delta_rt'].groupby(['peptide_vc']).describe()
    #for all non-variant findings
    nonvariant_stats=nonvariant_df['delta_rt'].groupby()
    #here plot interesting things
    overlap_stats.plot.hist()
    variant_df.groupby()
