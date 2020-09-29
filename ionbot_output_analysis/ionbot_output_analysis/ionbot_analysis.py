#!/usr/bin/env python3

import argparse, sys
import main_functions
import calculations
import helper_functions
import file_import
import plots
import fdr_reestimation
import pandas as pd
import logging
from Bio.SubsMat import MatrixInfo


'''
main analysis script.
'''

def main():
    '''this is the main function that will iterate over the giant pandas df and perform all analyses and make all figures
    this is written in a way that iteration should only be done once.
    '''
    parser = argparse.ArgumentParser(description='Ionbot output analysis')
    parser.add_argument('--ont', help='Directory ONT ionbot output files', required=True)
    parser.add_argument('--ref', help='Directory GENCODE reference ionbot output files', required=True)
    parser.add_argument('--cvc', help='Directory combi variant-containing', required=True)
    parser.add_argument('--cvf', help='Directory combi variant-free', required=True)
    parser.add_argument('--var', help='csv file of selected SAAV peptides', required=True)
    parser.add_argument('--decoy', help='Decoy peptide candidates for FDR re-estimation for variant-containing search',required=True)
    parser.add_argument('--bed', help='Bed file ONT isoforms', required=True)
    parser.add_argument('--gff', help='Gff3 file GENCODE isoforms', required=True)
    parser.add_argument('--rt', help='retention time (obs&pred)',required=True)
    parser.add_argument('--crap', help='Contamination file',required=True)
    args = vars(parser.parse_args()) 

    #import ionbot output data
    print("importing ionbot results")
    # overlap_prots = [line.strip() for line in open(args['ov'], 'r')]
    contam=helper_functions.prepare_contaminants(args['crap'])
    ibdf_ontonly=file_import.concatenate_csvs(args['ont'],contam,'ont')
    ibdf_refonly=file_import.concatenate_csvs(args['ref'],contam,'ref')
    ibdf_vf=file_import.concatenate_csvs(args['cvf'],contam,'vf')
    ibdf_vc=file_import.concatenate_csvs(args['cvc'],contam,'vc')
    
    #inital QC
    print("plotting initial QC")
    plots.plot_qvalues_comparison({'ONT only':ibdf_ontonly.dropna(subset=['q_value','DB']),'Ref only':ibdf_refonly.dropna(subset=['q_value','DB']),'Combi variant-containing':ibdf_vc.dropna(subset=['q_value','DB']),'Combi variant-free':ibdf_vf.dropna(subset=['q_value','DB'])},fdr_levels=[0.01])
    plots.plot_scores(ibdf_ontonly.dropna(),ibdf_refonly.dropna(),ibdf_vf.dropna())
    plots.plot_scores_combi(ibdf_vf.dropna(),ibdf_vc.dropna())
    plots.plot_target_decoy(ibdf_vf.dropna(),"qc_score_decoy_varfree.png", plot_title="Search result variant-free")
    plots.plot_target_decoy(ibdf_vc.dropna(),"qc_score_decoy_varcont.png", plot_title="Search result variant-containing")

    #read out paper-relevant information
    print('unexpected mofidication ratios in the VF set:')
    print(ibdf_vf['unexpected_modification'].value_counts(normalize=True))

    #import other data
    print('importing helper data')
    variant_peptides=file_import.il_sensitive_read_csv(args['var']).drop(columns=['id','haplotype']).drop_duplicates()
    decoy_variants=file_import.il_sensitive_read_csv(args['decoy']).drop(columns=['id']).drop_duplicates()
    rt_df=pd.read_csv(args['rt']).rename({'seq':'matched_peptide'},axis=1)
    
    #collect results
    print("Doing general analysis...")
    ###gather information about all non-variant matches###
    all_matches_nonvar_vf=ibdf_vf[(ibdf_vf['DB']==False)&(ibdf_vf['best_psm']==1)&(ibdf_vf['q_value']<0.01)].merge(rt_df[['matched_peptide','predicted_tr']],on='matched_peptide') #observed matches
    all_matches_nonvar_vc=ibdf_vc[(ibdf_vc['DB']==False)&(ibdf_vc['best_psm']==1)&(ibdf_vc['q_value']<0.01)].merge(rt_df[['matched_peptide','predicted_tr']],on='matched_peptide')
    all_match_unique_vf=all_matches_nonvar_vf.drop_duplicates(subset=['peptide'])
    all_match_unique_vc=all_matches_nonvar_vc.drop_duplicates(subset=['peptide'])
        
    # #get strand and chrom info
    print('...fetching and plotting origin info...')
    main_functions.origin_info_fetch(all_match_unique_vf,all_match_unique_vc,args['gff'],args['bed'])
    #get support
    print('...gauging protein support...')
    main_functions.protein_support(all_matches_nonvar_vf['proteins'],all_matches_nonvar_vc['proteins'])
    #get source
    print('...summing source dictionaries...')
    all_matches_nonvar_vc[all_matches_nonvar_vc['source_dict']=='ref'].drop_duplicates(subset='peptide').to_csv('reference_only_peptides.csv',index=False)
    main_functions.dict_source_bin(ibdf_vf['source_dict'],ibdf_vc['source_dict'],"sources_psm_theoretical_varfree.png","sources_psm_theoretical_varcont.png") #what we expect from the data
    main_functions.dict_source_bin(all_matches_nonvar_vf['source_dict'],all_matches_nonvar_vc['source_dict'],"sources_psm_obs_varfree.png","sources_psm_obs_varcont.png") #what we see
    main_functions.dict_source_bin(all_match_unique_vf['source_dict'],all_match_unique_vc['source_dict'],"sources_psm_obs_uniq_vf.png","sources_psm_obs_uniq_vc.png") #what we see
    # coverage - need to import full sequences for this #TODO
    # this is easier in the new ionbot output since the position numbers are included with the protein output
    # plots.plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varfree.png","vertical_coverage_varfree.png")

    #find variants
    print('Finding variants...')
    #for the variant-free set, recalculate FDR based on the subset of peptides that were predicted to have a SAAV
    detected_variant_combi_vf=ibdf_vf[(ibdf_vf["pred_aa_sub"]!="") & (ibdf_vf["DB"]==False)] #get all peptides that are predicted to have a SAAV
    detected_nonvar_combi_vf=ibdf_vf[(ibdf_vf["pred_aa_sub"]=="") & (ibdf_vf["DB"]==False)] #all not having saav
    detected_decoy_combi_vf=ibdf_vf[(ibdf_vf["pred_aa_sub"]!="") & (ibdf_vf["DB"]==True)]
    detected_nonvar_decoy_combi_vf=ibdf_vf[(ibdf_vf["pred_aa_sub"]=="") & (ibdf_vf["DB"]==True)]
    
    #for the variant-containing set, recalculate FDR based on the variant subset only
    observed_variants_vc=ibdf_vc.merge(variant_peptides.rename({'variant_peptide':'peptide'},axis=1), on='peptide').rename({'peptide':'variant_peptide'},axis=1) # target variants vc
    observed_decoy_vc=ibdf_vc.merge(decoy_variants.rename({'variant_peptide':'peptide'},axis=1), on='peptide').rename({'peptide':'variant_peptide'},axis=1) # decoy variants vc
    observed_decoy_vc=observed_decoy_vc[observed_decoy_vc["DB"]==True]
    observed_variant_counterparts_vc=ibdf_vc.merge(variant_peptides.rename({'ref_counterpart':'peptide'},axis=1), on='peptide').rename({'peptide':'ref_counterpart'},axis=1) # target counterparts vc
    observed_decoy_counterparts_vc=ibdf_vc.merge(decoy_variants.rename({'ref_counterpart':'peptide'},axis=1), on='peptide').rename({'peptide':'ref_counterpart'},axis=1) # decoy counterparts vc

    #just to get a count of how many variants were found in variant free set
    # print("variants found in variant-free before FDR correction:"+ )
    # print("variants found in the variant-containing before FDR correction:"+)
    print(str(detected_variant_combi_vf.drop_duplicates(subset='matched_peptide').shape[0])+' variants found in the variant-free output and '+str(observed_variants_vc.drop_duplicates(subset='matched_peptide').shape[0])+' variants found in the variant-containing output before FDR correction.')

    ###FDR###
    print('...filtering false positives...')
    #for variant-free set, filter for true variant peptides after the FDR
    prelim_variantset_vf=calculations.fdr_recalc_variantpep(detected_variant_combi_vf,detected_decoy_combi_vf,'variant_free_ppplot.png') # variants vf
    prelim_variantset_vf.merge(ibdf_vc,on='title',suffixes=('_vf','_vc')).to_csv('prelim_vf.csv',index=False)
    prelim_counterpartset_vf=calculations.fdr_recalc_variantpep(detected_nonvar_combi_vf,detected_nonvar_decoy_combi_vf,'variant_free_ctp_ppplot.png') # counterparts vf
    final_variantset_vf=prelim_variantset_vf.merge(variant_peptides.rename({'ref_counterpart':'peptide'},axis=1),how='left',on='peptide',indicator=True).rename({'peptide':'ref_counterpart'},axis=1).merge(rt_df,how='left',on='matched_peptide')
    print(f"VF variant peptides after FDR correction but before correction filter={prelim_variantset_vf.drop_duplicates().shape[0]}")
    final_variantset_vf=final_variantset_vf[final_variantset_vf['_merge']=='both'].drop(columns=['_merge']).reset_index(drop=True)
    print(f"out of all {str(final_variantset_vf.shape[0])} variant peptides correctly detected by ionbot, {final_variantset_vf[final_variantset_vf['substitution']==final_variantset_vf['pred_aa_sub']].shape[0]} of them had the correct aa substitution")
    final_counterpartset_vf=prelim_counterpartset_vf.merge(variant_peptides.rename({'ref_counterpart':'peptide'},axis=1)).rename({'peptide':'ref_counterpart'},axis=1).merge(rt_df,how='left',on='matched_peptide')#.merge(rt_pred_df,how='left',on='matched_peptide')

    #for the variant-free set, get true variant peptides from the FDR re-estimation
    final_variantset_vc=calculations.fdr_recalc_variantpep(observed_variants_vc,observed_decoy_vc,'variant_cont_ppplot.png').merge(rt_df,how='left',on='matched_peptide')
    final_counterpartset_vc=calculations.fdr_recalc_variantpep(observed_variant_counterparts_vc,observed_decoy_counterparts_vc,'variant_cont_ctp_ppplot.png').merge(rt_df,how='left',on='matched_peptide')
    print(str(final_variantset_vf.drop_duplicates(subset='matched_peptide').shape[0])+' variants found in the variant-free output and '+str(final_variantset_vc.drop_duplicates(subset='matched_peptide').shape[0])+' variants found in the variant-containing output after FDR correction.')
    final_variantset_vc.to_csv('variants_vc.csv',index=False)
    final_counterpartset_vc.to_csv('variants_vc_ctp.csv',index=False)
    final_variantset_vf.to_csv('variants_vf.csv',index=False)
    final_counterpartset_vf.to_csv('variants_vf_ctp.csv',index=False)

    ##find variant peptides in the lower-ranked 

    print("Making final plots...")
    main_functions.discrepancy_check(final_variantset_vf, final_variantset_vc,all_matches_nonvar_vf,all_matches_nonvar_vc,variant_peptides)
    # plots.plot_final_venns(final_variantset_vc,final_variantset_vf)
    return("Finished")
    
if __name__ == "__main__":
    main()
    
    