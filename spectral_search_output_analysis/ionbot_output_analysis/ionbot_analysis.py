#!/usr/bin/env python3

import argparse
import main_functions
import calculations
import helper_functions
import file_import
import plots


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

    #qc function
    print("plotting initial QC")
    plots.plot_scores(ibdf_ontonly.dropna(),ibdf_refonly.dropna(),ibdf_combi.dropna())
    plots.plot_scores_pg(ibdf_combi.dropna(),ibdf_combi_pg.dropna())
    plots.plot_scores_decoy(ibdf_combi.dropna(),"qc_pearsonr_decoy_varfree.png")
    plots.plot_scores_decoy(ibdf_combi_pg.dropna(),"qc_pearsonr_decoy_varcont.png")

    #filter badly scoring hits
    ibdf_ontonly = file_import.chunk_preprocessing(ibdf_ontonly)
    ibdf_refonly = file_import.chunk_preprocessing(ibdf_refonly)
    ibdf_combi = file_import.chunk_preprocessing(ibdf_combi)
    ibdf_combi_pg= file_import.chunk_preprocessing(ibdf_combi_pg)

    #import other data
    print('importing helper data')
    cpdt_pep,full_seqs=file_import.import_cpdt(args['cpdtvf'],True) #import cpdt will all non-variant-containing peptides (cat gencode and flair beforehand). full seqs for calculating horizontal coverage
    mut_cpdt,mut_pep_probs=file_import.import_cpdt(args['cpdtvar'],False) #import the cpdt file with all snv peptides
    mut_cpdt_counterparts,counterpart_pep_probs=file_import.import_cpdt(args['cpdtctp'],False)
    chromdict,stranddict=file_import.create_chromosome_reference(args['gff'],args['bed']) #import information about the chromosome of origin (QC)
    
    #iterate to fill the data structures
    print("doing analysis...")
    theoretical_saav= calculations.theoretical_saav_counts(mut_cpdt,mut_cpdt_counterparts)
    mut_observed_openmut,mutprotset,chromdist_openmut,stranddist_openmut=main_functions.combidict_analysis(ibdf_combi,chromdict,stranddict,cpdt_pep,theoretical_saav,mut_pep_probs,full_seqs,mut_cpdt,mut_cpdt_counterparts,True)
    mut_observed_classic,chromdist_classic,stranddist_classic=main_functions.combidict_analysis(ibdf_combi_pg,chromdict,stranddict,cpdt_pep,theoretical_saav,mut_pep_probs,full_seqs,mut_cpdt,mut_cpdt_counterparts,False)
    main_functions.discrepancy_check(mut_observed_classic,mut_observed_openmut, ibdf_combi, ibdf_combi_pg)
    plots.plot_chromosomal_dist(chromdist_classic,chromdist_openmut)
    plots.plot_strand_dist(stranddist_classic,stranddist_openmut)
    plots.plot_final_venns(mut_observed_classic,mut_observed_openmut,mut_cpdt,mutprotset)
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
args = vars(parser.parse_args()) 
main(args)
    
    
    
    