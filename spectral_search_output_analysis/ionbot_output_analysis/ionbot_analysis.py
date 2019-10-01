#!/usr/bin/env python3

import argparse
import main_functions
import calculations
import helper_functions
import file_import
import plots
import fdr_reestimation


'''
this script will concatenate all the separate spectral match files generated from ionbot and produce figures to answer research questions
'''

def main(args):
    '''this is the main function that will iterate over the giant pandas df and perform all analyses and make all figures
    this is written in a way that iteration should only be done once.
    '''
    #import ionbot output data
    print("importing ionbot results")
    ibdf_ontonly=file_import.concatenate_csvs(args['ont'])
    ibdf_refonly=file_import.concatenate_csvs(args['ref'])
    ibdf_combi=file_import.concatenate_csvs(args['cvf'])
    ibdf_combi_pg=file_import.concatenate_csvs(args['cvc'])

    #inital QC
    print("plotting initial QC")
    plots.plot_scores(ibdf_ontonly.dropna(),ibdf_refonly.dropna(),ibdf_combi.dropna())
    plots.plot_scores_pg(ibdf_combi.dropna(),ibdf_combi_pg.dropna())
    plots.plot_target_decoy(ibdf_combi.dropna(),"qc_pearsonr_decoy_varfree.png", plot_title="Search result variant-free")
    plots.plot_target_decoy(ibdf_combi_pg.dropna(),"qc_pearsonr_decoy_varcont.png", plot_title="Search result variant-containing")
    plots.plot_qvalues_comparison({'ONT only':ibdf_ontonly,'Ref only':ibdf_refonly,'Combi variant-containing':ibdf_combi_pg,'Combi variant-free':ibdf_combi},fdr_levels=[0.01])


    #filter badly scoring hits and decoys - now incorporating this filtering in the main analysis
    # ibdf_ontonly = file_import.post_filtering(ibdf_ontonly_original)
    # ibdf_refonly = file_import.post_filtering(ibdf_refonly_original)
    # ibdf_combi = file_import.post_filtering(ibdf_combi_original)
    # ibdf_combi_pg= file_import.post_filtering(ibdf_combi_pg_original)

    ##########
    #idea: move data into SQL database and move the quality control plots above to another script to cut down on processing time, since the quality control steps are the only ones that need all, unfiltered data
    #insert here an SQL query here to import data that meets preprocessing cutoffs
    #########

    #import other data
    print('importing helper data')
    cpdt_pep,full_seqs=file_import.import_cpdt(args['cpdtvf'],True) #import cpdt will all non-variant-containing peptides (cat gencode and flair beforehand). full seqs for calculating horizontal coverage
    mut_cpdt,mut_pep_probs=file_import.import_cpdt(args['cpdtvar'],False) #import the cpdt file with all snv peptides
    mut_cpdt_counterparts,counterpart_pep_probs=file_import.import_cpdt(args['cpdtctp'],False)
    chromdict,stranddict=file_import.create_chromosome_reference(args['gff'],args['bed']) #import information about the chromosome of origin (QC)
    cpdt_rev=file_import.import_cpdt_simple(args['decoy']) #import information about variant decoys
    cpdt_rev_ctp=file_import.import_cpdt_simple(args['decoyctp'])
    # rt=pd.read_csv(args['rt'])
    
    #iterate to fill the data structures
    print("doing analysis...")
    theoretical_saav= calculations.theoretical_saav_counts(mut_cpdt,mut_cpdt_counterparts)
    mut_observed_openmut,mutprotset,chromdist_openmut,stranddist_openmut=main_functions.combidict_analysis(ibdf_combi,chromdict,stranddict,cpdt_pep,full_seqs,theoretical_saav,mut_pep_probs,mut_cpdt,mut_cpdt_counterparts,cpdt_rev,cpdt_rev_ctp,True)
    mut_observed_classic,chromdist_classic,stranddist_classic=main_functions.combidict_analysis(ibdf_combi_pg,chromdict,stranddict,cpdt_pep,full_seqs,theoretical_saav,mut_pep_probs,mut_cpdt,mut_cpdt_counterparts,cpdt_rev,cpdt_rev_ctp,False)
    #here do a re-calculation of FDR for variant peps (look into what functions from combidict_analysis need to be moved to after this recalculation)
    # print("recalculating FDR for variant peptides...")
    # print(mut_observed_classic)
    # print(str(sum(mut_observed_classic.values()))+' variant containing and '+str(sum(mut_observed_openmut.values()))+' variant free observed variants before FDR correction...')
    # var_vf,var_vc=fdr_reestimation.main(ibdf_combi.dropna(),ibdf_combi_pg.dropna(),raw_vc,mut_observed_classic,args["decoy"])
    
    #get variant groupings
    abrv_vc,ref_abrv_vc=helper_functions.abbreviate_peps(mut_observed_openmut)
    abrv_vf,ref_abrv_vf=helper_functions.abbreviate_peps(mut_observed_classic)
    print("making final plots...")
    main_functions.discrepancy_check(ref_abrv_vc,ref_abrv_vf,var_vc,var_vf, ibdf_combi, ibdf_combi_pg,args['rt'])
    plots.plot_chromosomal_dist(chromdist_classic,chromdist_openmut)
    plots.plot_strand_dist(stranddist_classic,stranddist_openmut)
    plots.plot_final_venns(abrv_vc,abrv_vf,mut_cpdt,mutprotset)
    return("Finished")
    
parser = argparse.ArgumentParser(description='Ionbot output analysis')
parser.add_argument('--ont', help='Directory ONT ionbot output files', required=True)
parser.add_argument('--ref', help='Directory GENCODE reference ionbot output files', required=True)
parser.add_argument('--cvc', help='Directory combi variant-containing', required=True)
parser.add_argument('--cvf', help='Directory combi variant-free', required=True)
parser.add_argument('--cpdtvar', help='CPDT file of selected SAAV peptides', required=True)
parser.add_argument('--cpdtctp', help='CPDT file of non-mutated counterparts of the selected SAAV peptides', required=True)
parser.add_argument('--cpdtvf', help='CPDT file of entire combi variant-free', required=True)
parser.add_argument('--bed', help='Bed file ONT isoforms', required=True)
parser.add_argument('--gff', help='Gff3 file GENCODE isoforms', required=True)
parser.add_argument('--decoy', help='Decoy peptide candidates for FDR re-estimation for variant-containing search',required=True)
parser.add_argument('--decoyctp', help='Decoy counterpart peptide candidates for FDR re-estimation for variant-containing search',required=True)
parser.add_argument('--rt', help='retention time prediction',required=True)
# parser.add_argument('--varpeps' help='Scan IDs of variant peptides that were identified as "true" variant peptides')
args = vars(parser.parse_args()) 
main(args)
    
    
    
    