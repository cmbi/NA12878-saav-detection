#!/usr/bin/env python3

import re, os,sys, itertools, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
import plots
import helper_functions
import calculations

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
    origin_info=file_import.create_chromosome_reference(gff,bed)
    vf=df_vf.merge(origin_info,on='transcript_id') #this came from taking the first protein on the list and getting its transcript id
    vc=df_vc.merge(origin_info,on='transcript_id')
    strand_vf=Counter(dict(vf['strand'].value_counts()))
    strand_vc=Counter(dict(vc['strand'].value_counts()))
    chrom_vf=Counter(dict(vf['chromosome'].value_counts()))
    chrom_vc=Counter(dict(vc['chromosome'].value_counts()))
    plots.plot_chromosomal_dist(chrom_vc,chrom_vf)
    plots.plot_strand_dist(strand_vc,strand_vf)

def dict_source_bin(df_vf,df_vc,figname_vf,figname_vc):
    vf=Counter(dict(df_vf.value_counts()))
    vc=Counter(dict(df_vc.value_counts()))
    plots.plot_source_piechart(vf,figname_vf)
    plots.plot_source_piechart(vc,figname_vc)        

def discrepancy_check(vf,vc,nonvar_vf,nonvar_vc,rt):
    '''check out the differences in identifications between the 2 combination dictionaries
    why doesn't ionbot catch everything? look at the ones that it does not catch but the variant-containing dictionary does
    plot lengths of the missed/caught peptides (longer than average?)
    plot unexpected modifications of the missed peptides (more unexpected modifications than average?)'''
    #abbreviate the results, group by matched peptide
    aggregation_rules={'percolator_psm_score':'mean','unexpected_modification':helper_functions.longest,'proteins':helper_functions.longest,'matched_peptide':helper_functions.longest}
    ibdf_vf=vf.groupby('peptide').aggregate(aggregation_rules)
    ibdf_vc=vc.groupby('peptide').aggregate(aggregation_rules)
    #convert everything to sets and find the variants that are/are not in common
    merged=pd.merge(ibdf_vf,ibdf_vc, on='peptide', how='outer', suffixes=('_vf','_vc'), indicator=True)
    vf_only=merged[merged['_merge'] == 'left_only']
    vc_only=merged[merged['_merge'] == 'right_only']
    overlap=merged[merged['_merge'] == 'both']
    #plot lengths of the peptides caught by one method but not the other
    plots.plot_peplengths(Counter(vc_only['matched_peptide_vc'].str.len().to_list()),Counter(vf_only['matched_peptide_vf'].str.len().to_list()),variant=True)
    plots.plot_peplengths(Counter(nonvar_vc['matched_peptide'].str.len().to_list()),Counter(nonvar_vf['matched_peptide'].str.len().to_list()),variant=False)
    #look at unexpected modifications in the mis-labeled VF (since there are no VF exclusive variants)
    vc_unique=pd.merge(vf['peptide'],vc[['scan_id','peptide','matched_peptide','percolator_psm_score']], on='peptide', how='right', indicator=True) #isolate variant containing only
    mislabeled=pd.merge(nonvar_vf[['scan_id','peptide','unexpected_modification']],vc_unique.loc[vc_unique['_merge']=='right_only','scan_id'], on='scan_id').groupby(['peptide','unexpected_modification']).aggregate(lambda x: x.iloc[0]).reset_index()
    plots.plot_unexpected_mods(mislabeled["unexpected_modification"].apply(helper_functions.categorize_mods),pd.DataFrame(),variant=True) #plot what modifications they contained
    plots.plot_unexpected_mods(nonvar_vf["unexpected_modification"].apply(helper_functions.categorize_mods),nonvar_vc["unexpected_modification"].apply(helper_functions.categorize_mods),variant=False) #plot what modifications they contained
    #direct comparison of scores of peptides
    mislabeled=pd.merge(nonvar_vf[['scan_id','matched_peptide','percolator_psm_score']],vc_unique.loc[vc_unique['_merge']=='right_only',('scan_id','matched_peptide','percolator_psm_score')], on='scan_id',suffixes=('_vf','_vc')).groupby(['matched_peptide_vc','matched_peptide_vf']).aggregate('mean').reset_index()
    scores_all_filtered=mislabeled.loc[mislabeled["percolator_psm_score_vc"]>mislabeled["percolator_psm_score_vf"]] #information for scan ids that are higher in variant containing than variant free
    scores_all_filtered.to_csv('vc_higher_than_vf.csv',index=False) #first print to csv, look at where the vc scores were higher than vf
    plots.plot_ib_scores_directcomp(mislabeled.rename(columns={"matched_peptide_vc":"matched_peptide"}),rt) #direct comparison plot: what scores they had in each of the libraries
    #general comparison of the scores
    plots.plot_ib_scores(vf_only["percolator_psm_score_vf"].tolist(),vc_only["percolator_psm_score_vc"].tolist(),overlap["percolator_psm_score_vc"].tolist(),overlap["percolator_psm_score_vf"].tolist(),nonvar_vc["percolator_psm_score"].tolist(),nonvar_vf["percolator_psm_score"].tolist())
    #to explore: return scan ids and check the ids that were not identified with variant free method in a later function
   