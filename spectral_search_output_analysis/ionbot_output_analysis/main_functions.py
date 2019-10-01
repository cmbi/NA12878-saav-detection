#!/usr/bin/env python3

import re, os,sys, itertools, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
import plots
import helper_functions
import calculations
 
def combidict_analysis(results_df,chromdict,stranddict,cpdt_pep,full_seqs,theoretical_saavs,mut_pep_probs,mut_cpdt_theoretical,mut_cpdt_counterparts,cpdt_decoyvar,cpdt_decoyctp,isOpenmut):
    ''' 
    the main analysis function of the package
    go through the results list (after minimal filtering) and add the relevant information to various objects


    '''
    # proteins_covered=Counter() #proteins detected
    mutated=set() #all proteins that were detected to have a variant by ionbot. how does compare to the proteins that actually do have variant?
    ref_only=set() #scan ids matching to proteins in the reference set
    ont_only=set() #scan ids matching to proteins in the ont set
    both=set() #scan ids that matched to both ref and ont proteins
    chrom_dist=Counter() # 1 chromosome location per scan id
    strand_dist=Counter() # associated strand
    target_frame=pd.DataFrame() # true variant peptides
    decoy_frame=pd.DataFrame() # decoy variant peptides (for purposes of FDR)
    counterpart_frame=pd.DataFrame() #true reference counterpart peptides
    decoy_ctp_frame=pd.DataFrame() # decoy reference counterpart peptides (for purposes of FDR)
    protein_support=Counter() # how many peptides per protein found
    unamb_protsupport=Counter() # how many unambiguous peptides per protein
    for row in tqdm(results_df.iterrows()): #iterate through the results dataframe
        scanid=row[1][0]
        mod=str(row[1][7])
        aamod=re.findall('[A-Z]->[A-Z]',mod)
        pep=row[1][3]
        prot_ids=row[1][9]
        q_val=row[1][11]
        decoy=row[1][6]
        if q_val<0.01:
            if '||' in prot_ids:
                ids=prot_ids.split('||')
                if not decoy:
                    for i in ids:
                        protein_support[helper_functions.get_id(i)]+=1
            else: #unambiguous assignment!
                ids=[prot_ids]
                if not decoy:
                    unamb_protsupport[helper_functions.get_id(prot_ids)]+=1
            if not decoy:
                # proteins_covered=detected_proteins(ids,proteins_covered) #what proteins from the proteome are covered and in what amounts
                cpdt_pep=helper_functions.add_to_observed(pep,ids,cpdt_pep,False) #what peptides are detected, how many, and what proteins they come from
                cptfound=helper_functions.detect_peptides(pep,ids,mut_cpdt_counterparts,isOpenmut) #note the occurences of "reference" versions of the SAAV peptides
                if cptfound:
                    counterpart_frame=counterpart_frame.append(results_df.loc[[row[0]]])
                chrom_origin=helper_functions.find_chrom(ids,chromdict)
                strand_origin=helper_functions.find_strand(ids,stranddict)
                chrom_dist[chrom_origin]+=1 #which chromosome does the peptide belong to
                strand_dist[strand_origin]+=1 #which strand does the peptide belong to
                ref_only,ont_only,both=helper_functions.bin_hits_by_source(scanid,ids,ref_only,ont_only,both, isOpenmut)#what dictionary/ies did the peptide match to
                if isOpenmut and len(aamod)>0: #if ib detects a mutated peptide (for open variant search only)
                    # mut_cpdt_observed,varfound=helper_functions.add_to_observed(pep,ids,mut_cpdt_observed,isOpenmut)
                    varfound=helper_functions.detect_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)
                    if varfound:
                        # observed_raw.add(pep)
                        target_frame=target_frame.append(results_df.loc[[row[0]]])
                        print(target_frame)
                        sys.exit()
                    for i in ids: 
                        mutated.add(helper_functions.get_id(i))
                # elif isOpenmut and detect_mut_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)!='': ##very strange scenario here!!##
                #     print(scanid)
                elif not isOpenmut: #check if mutant peptide if not open mutation settings
                    # mut_cpdt_observed,varfound=helper_functions.add_to_observed(pep,ids,mut_cpdt_observed,isOpenmut)
                    varfound=helper_functions.detect_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)
                    if varfound:
                        # observed_raw.add(pep)
                        target_frame=target_frame.append(results_df.loc[[row[0]]])
            elif decoy:
                cpdt=[]
                cpdt_ctp=[]
                for i in ids:
                    if i in cpdt_decoyvar:
                        cpdt+=cpdt_decoyvar[i]
                    elif i in cpdt_decoyctp:
                        cpdt_ctp+=cpdt_decoyctp[i]
                for v in cpdt:
                    if pep in v or helper_functions.equivalent_check(pep,v):
                        decoy_frame=decoy_frame.append(df_in.loc[[row[0]]])
                for c in cpdt_ctp:
                    if pep in c or helper_functions.equivalent_check(pep,c):
                        decoy_ctp_frame=decoy_ctp_frame.append(df_in.loc[[row[0]]])
    #do FDR re-estimation/filtering
    variant_frame=fdr_recalc_variantpep(target_frame,decoy_frame)
    counterpart_frame=fdr_recalc_variantpep(counterpart_frame,decoy_ctp_frame)
    #compile all the variants after filtering
    var_vf=Counter(dict(variant_frame['matched_peptide'].value_counts()))
    print(str(sum(var_vf.values()))+' variants remaining after FDR correction.')
    mut_cpdt_observed=df_to_dict(variant_frame,mut_cpdt_theoretical,isOpenmut)
    mut_counterparts_observed=df_to_dict(counterpart_frame,mut_cpdt_counterparts,isOpenmut)
    #some plots
    if isOpenmut:
        # with open('checkpoint.openmut.json','w'):
        #     json.dump()
        # plots.plot_mut_abundance(mut_cpdt_observed,mut_counterparts_observed,cpdt_pep,full_seqs,"mutant_abundance_varfree.png")
        plots.plot_mut_vs_nonmut(mut_cpdt_observed,mut_counterparts_observed,theoretical_saavs,mut_pep_probs,"_varfree.png")
        plots.plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varfree.png","vertical_coverage_varfree.png")
        plots.plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varfree.png",isOpenmut) #must stay in
        plots.plot_support(protein_support,unamb_protsupport,'protein_evidence_varfree.png') #must stay in, do we only want to look at support for variant peptides or keep global measures here?
    else:
        # plots.plot_mut_abundance(mut_cpdt_observed,mut_counterparts_observed,cpdt_pep,full_seqs,"mutant_abundance_varcont.png")
        plots.plot_mut_vs_nonmut(mut_cpdt_observed,mut_counterparts_observed,theoretical_saavs,mut_pep_probs,"_varcont.png")
        plots.plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varcont.png","vertical_coverage_varcont.png")
        plots.plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varcont.png",isOpenmut)
        plots.plot_support(protein_support,unamb_protsupport,'protein_evidence_varcont.png')
    if isOpenmut:
        return(helper_functions.count_muts(mut_cpdt_observed),mutated,chrom_dist,strand_dist)
    return(helper_functions.count_muts(mut_cpdt_observed),chrom_dist,strand_dist)

def fdr_recalc_variantpep(target,decoy):
    '''filter df by new q-value
    '''
    print('Recalculating q-values for variant peptide set...')
    print(str(target.shape[0])+' target variant peptides and '+str(decoy.shape[0])+' decoy variant peptides found.')
    df_in=pd.concat([target,decoy],ignore_index=True)
    df_in.columns=['scan_id','charge','precursor_mass','matched_peptide','modifications','percolator_psm_score','DB','unexpected_modification','ms2pip_pearsonr','proteins','num_unique_pep_ids','q_value']
    df=calculations.calculate_qvalues(df_in,decoy_col='DB',score_col='percolator_psm_score')
    indices=np.argwhere(df['q_value']<0.01)
    return(df[:int(indices[-1][0])+1]) #threshold cut

def df_to_dict(filtered_variants,mut_cpdt_theoretical,isOpenmut):
    ''' process the variant information post FDR re-estimation
    '''
    mut_cpdt_observed=mut_cpdt_theoretical
    variant_counterpart_couples=[]
    for row in filtered_variants.iterrows():
        pep=row[1][3]
        prot_ids=row[1][9]
        decoy=row[1][6]
        if '||' in prot_ids:
            ids=prot_ids.split('||')
        else:
            ids=[prot_ids]
        if not decoy:
            mut_cpdt_observed,varfound=helper_functions.add_to_observed(pep,ids,mut_cpdt_observed,isOpenmut)
            # mut_counterparts_observed,cptfound=helper_functions.add_to_observed(pep,ids,mut_counterparts_observed,isOpenmut,variant_check=True)
    return(helper_functions.remove_empty(mut_observed))
        

def discrepancy_check(dict_saavs_vc,dict_saavs_vf,allmuts_classic,allmuts_openmut,ibdf_combi,ibdf_combi_pg,rt):
    '''check out the differences in identifications between the 2 combination dictionaries
    why doesn't ionbot catch everything? look at the ones that it does not catch but the variant-containing dictionary does
    plot lengths of the missed/caught peptides (longer than average?)
    plot unexpected modifications of the missed peptides (more unexpected modifications than average?)'''
    #convert everything to sets and find the variants that are/are not in common
    combined_ref=helper_functions.concat_dicts(dict_saavs_vc,dict_saavs_vf)
    discrepancy=helper_functions.unpack_grouped_variants(set(dict_saavs_vc.keys()).difference(set(dict_saavs_vf.keys())),dict_saavs_vc)
    agreement=helper_functions.unpack_grouped_variants(set(dict_saavs_vc.keys()).intersection(set(dict_saavs_vf.keys())),combined_ref)
    union=helper_functions.unpack_grouped_variants(set(dict_saavs_vc.keys()).union(set(dict_saavs_vf.keys())),combined_ref)
    ibonly=helper_functions.unpack_grouped_variants(set(dict_saavs_vf.keys()).difference(set(dict_saavs_vc.keys())),dict_saavs_vf)
    discrepancy_ct_pg=allmuts_classic-allmuts_openmut
    discrepancy_ct_om=allmuts_openmut-allmuts_classic
    ibdf_combi_nonmut=ibdf_combi.loc[~ibdf_combi["matched_peptide"].isin(union)]
    ibdf_combi_pg_nonmut=ibdf_combi_pg.loc[~ibdf_combi_pg["matched_peptide"].isin(union)]
    #plot lengths of the peptides caught by one method but not the other
    non_mut_peplen_vf=Counter(dict(ibdf_combi_nonmut['matched_peptide'].value_counts()))
    non_mut_peplen_vc=Counter(dict(ibdf_combi_pg_nonmut['matched_peptide'].value_counts()))
    plots.plot_peplengths(discrepancy_ct_pg,discrepancy_ct_om,variant=True)
    plots.plot_peplengths(non_mut_peplen_vc,non_mut_peplen_vf,variant=False)
    #get the scan ids corresponding to all the variant peptides that are in one set but not the other
    scanids=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),"scan_id"].tolist()
    #Compare unexpected modifications
    unexp_mod_om=ibdf_combi.loc[ibdf_combi["scan_id"].isin(scanids),"unexpected_modification"].tolist()
    unexp_mod_pg=ibdf_combi_pg.loc[ibdf_combi_pg["scan_id"].isin(scanids),"unexpected_modification"].tolist()
    unexp_mod_nonmut_om=ibdf_combi_nonmut["unexpected_modification"].tolist()
    unexp_mod_nonmut_pg=ibdf_combi_pg_nonmut["unexpected_modification"].tolist()
    plots.plot_unexpected_mods(unexp_mod_om,unexp_mod_pg,variant=True) #plot what modifications they contained
    plots.plot_unexpected_mods(unexp_mod_nonmut_om,unexp_mod_nonmut_pg,variant=False) #plot what modifications they contained
    #direct comparison of scores of peptides found with variant-containing but not variant-free
    scores_pg=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),["scan_id","percolator_psm_score","matched_peptide"]]#.tolist()
    scores_om=ibdf_combi.loc[ibdf_combi["scan_id"].isin(scanids),["scan_id","percolator_psm_score"]] #added matched peptide for length dimension
    scores_all=scores_pg.merge(scores_om, on=('scan_id'),suffixes=('_vc','_vf'))
    scores_all=scores_all.loc[scores_all["percolator_psm_score_vc"]>scores_all["percolator_psm_score_vf"]] #information for scan ids that are higher in variant containing than variant free
    scores_all.to_csv('vc_higher_than_vf.csv',index=False) #first print to csv
    plots.plot_ib_scores_directcomp(scores_om,scores_pg,rt) #direct comparison plot: what scores they had in each of the libraries
    #general comparison of the scores
    list_ibonly=ibdf_combi.loc[ibdf_combi["matched_peptide"].isin(ibonly),"percolator_psm_score"].tolist()
    list_intersection_varfree=ibdf_combi.loc[ibdf_combi["matched_peptide"].isin(agreement),"percolator_psm_score"].tolist()
    list_nonmut_vf=ibdf_combi_nonmut["percolator_psm_score"].tolist()
    print(len(list_intersection_varfree))
    list_intersection_varcont=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(agreement),"percolator_psm_score"].tolist()
    list_nonmut_vc=ibdf_combi_pg_nonmut["percolator_psm_score"].tolist()
    print(len(list_intersection_varcont))
    list_pgonly=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),"percolator_psm_score"].tolist()
    plots.plot_ib_scores(list_ibonly,list_pgonly,list_intersection_varcont,list_intersection_varfree,list_nonmut_vc,list_nonmut_vf)
    #to explore: return scan ids and check the ids that were not identified with variant free method in a later function
    return(0)
   