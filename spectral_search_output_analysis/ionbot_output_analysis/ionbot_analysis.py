#!/usr/bin/env python3

import argparse
import main_functions
import calculations
import helper_functions
import file_import
import plots
import fdr_reestimation
import pandas as pd


'''
this script will concatenate all the separate spectral match files generated from ionbot and produce figures to answer research questions
'''

def main(args):
    '''this is the main function that will iterate over the giant pandas df and perform all analyses and make all figures
    this is written in a way that iteration should only be done once.
    '''
    #import ionbot output data
    print("importing ionbot results")
    # ibdf_ontonly=file_import.concatenate_csvs(args['ont'])
    # ibdf_refonly=file_import.concatenate_csvs(args['ref'])
    ibdf_vf=file_import.concatenate_csvs(args['cvf'])
    ibdf_vc=file_import.concatenate_csvs(args['cvc'])

    #inital QC
    print("plotting initial QC")
    # plots.plot_scores(ibdf_ontonly.dropna(),ibdf_refonly.dropna(),ibdf_vf.dropna())
    plots.plot_scores_combi(ibdf_vf.dropna(),ibdf_vc.dropna())
    plots.plot_target_decoy(ibdf_vf.dropna(),"qc_pearsonr_decoy_varfree.png", plot_title="Search result variant-free")
    plots.plot_target_decoy(ibdf_vc.dropna(),"qc_pearsonr_decoy_varcont.png", plot_title="Search result variant-containing")
    # plots.plot_qvalues_comparison({'ONT only':ibdf_ontonly,'Ref only':ibdf_refonly,'Combi variant-containing':ibdf_vc,'Combi variant-free':ibdf_vf},fdr_levels=[0.01])

    #import other data
    print('importing helper data')
    origin_info=file_import.create_chromosome_reference(args['gff'],args['bed']) #import information about the chromosome of origin (QC), dataframe
    all_peptides=file_import.il_sensitive_read_csv(args['vfd'])
    variant_peptides=file_import.il_sensitive_read_csv(args['var'])
    variant_counterparts=file_import.il_sensitive_read_csv(args['ctp'])
    decoy_variants=file_import.il_sensitive_read_csv(args['decoy'])
    decoy_counterparts=file_import.il_sensitive_read_csv(args['decoyctp'])
    theoretical_saav= calculations.saav_counts(variant_peptides,variant_counterparts)
    
    #collect results
    print("Doing general analysis...")
    ###gather information about all non-variant matches###
    all_matches_nonvar_vf=ibdf_vf.merge(all_peptides, on='peptide') #observed matches to variant free dictionary
    all_matches_nonvar_vc=ibdf_vc.merge(all_peptides, on='peptide') #observed matches to variant containing dictionary
    #get strand and chrom info
    print('...fetching and plotting origin info...')
    main_functions.origin_info_fetch(all_matches_nonvar_vf,all_matches_nonvar_vc,origin_info)
    #get support
    print('...gauging protein support...')
    main_functions.protein_support(all_matches_nonvar_vf['proteins'],all_matches_nonvar_vc['proteins'])
    #get source
    print('...summing source dictionaries...')
    main_functions.dict_source_bin(ibdf_vf['source_dict'],ibdf_vc['source_dict'],"sources_psm_theoretical_varfree.png","sources_psm_theoretical_varcont.png") #what we expect from the data
    main_functions.dict_source_bin(all_matches_nonvar_vf['source_dict'],all_matches_nonvar_vc['source_dict'],"sources_psm_obs_varfree.png","sources_psm_obs_varcont.png") #what we see
    #coverage - need to import full sequences for this #TODO
    # plots.plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varfree.png","vertical_coverage_varfree.png")

    #find variants
    print('Finding variants...')
    #for the variant-free set, recalculate FDR based on the subset of peptides that were predicted to have a SAAV
    detected_variant_combi_vf=ibdf_vf[(ibdf_vf["unexpected_modification"].str.contains('[A-Z]->[A-Z]',regex=True)) & (ibdf_vf["DB"]==False)] #get all peptides that are predicted to have a SAAV
    detected_normal_combi_vf=ibdf_vf[(~ibdf_vf["unexpected_modification"].str.contains('[A-Z]->[A-Z]',regex=True)) & (ibdf_vf["DB"]==False)] #all not having saav
    detected_decoy_combi_vf=ibdf_vf[(ibdf_vf["unexpected_modification"].str.contains('[A-Z]->[A-Z]',regex=True)) & (ibdf_vf["DB"]==True)]
    detected_normal_decoy_combi_vf=ibdf_vf[(~ibdf_vf["unexpected_modification"].str.contains('[A-Z]->[A-Z]',regex=True)) & (ibdf_vf["DB"]==True)]

    # observed_variants_vf=ibdf_vf.merge(variant_peptides, ) # target variants vf
    # observed_decoy_vf=ibdf_vf.merge(decoy_variants, on='peptide') # decoy variants vf
    # observed_decoy_vf=observed_decoy_vf[observed_decoy_vf["DB"]==True]
    # observed_variant_counterparts_vf=ibdf_vf.merge(variant_counterparts, on='peptide') # target counterparts vf
    # observed_decoy_counterparts_vf=ibdf_vf.merge(decoy_counterparts, on='peptide') # decoy counterparts vf
    
    #for the variant-containing set, recalculate FDR based on the variant subset only
    observed_variants_vc=ibdf_vc.merge(variant_peptides, on='peptide') # target variants vc
    observed_decoy_vc=ibdf_vc.merge(decoy_variants, on='peptide') # decoy variants vc
    observed_decoy_vc=observed_decoy_vc[observed_decoy_vc["DB"]==True]
    observed_variant_counterparts_vc=ibdf_vc.merge(variant_counterparts, on='peptide') # target counterparts vc
    observed_decoy_counterparts_vc=ibdf_vc.merge(decoy_counterparts, on='peptide') # decoy counterparts vc

    ###FDR###
    print('...filtering false positives...')
    #for variant-free set, filter for true variant peptides after the FDR
    prelim_variantset_vf=calculations.fdr_recalc_variantpep(detected_variant_combi_vf,detected_decoy_combi_vf,'variant_free_ppplot.png') # variants vf
    prelim_counterpartset_vf=calculations.fdr_recalc_variantpep(detected_normal_combi_vf,detected_normal_decoy_combi_vf,'variant_free_ctp_ppplot.png') # counterparts vf
    final_variantset_vf=prelim_variantset_vf.merge(variant_peptides,on='peptide')
    final_counterpartset_vf=prelim_counterpartset_vf.merge(variant_counterparts,on='peptide')

    #for the variant-free set, get true variant peptides from the FDR re-estimation
    final_variantset_vc=calculations.fdr_recalc_variantpep(observed_variants_vc,observed_decoy_vc,'variant_cont_ppplot.png')
    final_counterpartset_vc=calculations.fdr_recalc_variantpep(observed_variant_counterparts_vc,observed_decoy_counterparts_vc,'variant_cont_ctp_ppplot.png')
    print(str(final_variantset_vf.shape[0])+' variants found in the variant-free output and '+str(final_variantset_vc.shape[0])+' variants found in the variant-containing output.')

    print('Analyzing variants...')
    sub_type_vc,sub_count_vc=calculations.saav_counts(final_variantset_vc,final_counterpartset_vc,observed=True)
    sub_type_vf,sub_count_vf=calculations.saav_counts(final_variantset_vf,final_counterpartset_vf,observed=True)
    plots.plot_heatmaps(theoretical_saav,'heatmap_theoretical_subs.png')
    plots.plot_heatmaps(sub_type_vc,'heatmap_obs_subs_vc.png')
    plots.plot_heatmaps(sub_type_vf,'heatmap_obs_subs_vf.png')
    plots.plot_mut_vs_nonmut(sub_count_vc,'variant_vs_counterpart_vc.png')
    plots.plot_mut_vs_nonmut(sub_count_vf,'variant_vs_counterpart_vf.png')

    print("Making final plots...")
    main_functions.discrepancy_check(final_variantset_vf, final_variantset_vc,all_matches_nonvar_vf,all_matches_nonvar_vc, args['rt'])
    plots.plot_final_venns(final_variantset_vc,final_variantset_vf)
    return("Finished")
    
parser = argparse.ArgumentParser(description='Ionbot output analysis')
parser.add_argument('--ont', help='Directory ONT ionbot output files', required=True)
parser.add_argument('--ref', help='Directory GENCODE reference ionbot output files', required=True)
parser.add_argument('--cvc', help='Directory combi variant-containing', required=True)
parser.add_argument('--cvf', help='Directory combi variant-free', required=True)
parser.add_argument('--var', help='csv file of selected SAAV peptides', required=True)
parser.add_argument('--ctp', help='csv file of non-mutated counterparts of the selected SAAV peptides', required=True)
parser.add_argument('--vfd', help='csv file of entire combi variant-free', required=True)
parser.add_argument('--bed', help='Bed file ONT isoforms', required=True)
parser.add_argument('--gff', help='Gff3 file GENCODE isoforms', required=True)
parser.add_argument('--decoy', help='Decoy peptide candidates for FDR re-estimation for variant-containing search',required=True)
parser.add_argument('--decoyctp', help='Decoy counterpart peptide candidates for FDR re-estimation for variant-containing search',required=True)
parser.add_argument('--rt', help='retention time prediction',required=True)
# parser.add_argument('--varpeps' help='Scan IDs of variant peptides that were identified as "true" variant peptides')
args = vars(parser.parse_args()) 
main(args)
    
    
    
    