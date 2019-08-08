#!/usr/bin/env python3

import re, os,sys, itertools, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
import plots
import helper_functions


def discrepancy_check(allmuts_classic,allmuts_openmut,ibdf_combi,ibdf_combi_pg):
    '''check out the differences in identifications between the 2 combination dictionaries
    why doesn't ionbot catch everything? look at the ones that it does not catch but the variant-containing dictionary does
    plot lengths of the missed/caught peptides (longer than average?)
    plot unexpected modifications of the missed peptides (more unexpected modifications than average?)'''
    discrepancy=set(allmuts_classic).difference(set(allmuts_openmut))
    agreement=set(allmuts_classic).intersection(set(allmuts_openmut))
    ibonly=set(allmuts_openmut).difference(set(allmuts_classic))
    discrepancy_ct_pg=allmuts_classic-allmuts_openmut
    discrepancy_ct_om=allmuts_openmut-allmuts_classic
    #plot lengths of the peptides caught by one method but not the other
    plots.plot_peplengths(discrepancy_ct_pg,discrepancy_ct_om)
    #get the scan ids corresponding to all the variant peptides that are in one set but not the other
    scanids=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),"scan_id"].tolist()
    #Compare unexpected modifications
    unexp_mod_om=ibdf_combi.loc[ibdf_combi["scan_id"].isin(scanids),"unexpected_modification"].tolist()
    unexp_mod_pg=ibdf_combi_pg.loc[ibdf_combi_pg["scan_id"].isin(scanids),"unexpected_modification"].tolist()
    plots.plot_unexpected_mods(unexp_mod_om,unexp_mod_pg) #plot what modifications they contained
    #direct comparison of scores of peptides found with variant-containing but not variant-free
    scores_pg=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),["scan_id","percolator_psm_score","matched_peptide"]]#.tolist()
    scores_om=ibdf_combi.loc[ibdf_combi["scan_id"].isin(scanids),["scan_id","percolator_psm_score"]] #added matched peptide for length dimension
    scores_all=scores_pg.merge(scores_om, on=('scan_id'),suffixes=('_vc','_vf'))
    scores_all=scores_all.loc[scores_all["percolator_psm_score_vc"]>scores_all["percolator_psm_score_vf"]] #information for scan ids that are higher in variant containing than variant free
    scores_all.to_csv('vc_higher_than_vf.csv',index=False) #first print to csv
    plots.plot_ib_scores_directcomp(scores_om,scores_pg) #direct comparison plot: what scores they had in each of the libraries
    #general comparison of the scores
    list_ibonly=ibdf_combi.loc[ibdf_combi["matched_peptide"].isin(ibonly),"percolator_psm_score"].tolist()
    list_intersection_varfree=ibdf_combi.loc[ibdf_combi["matched_peptide"].isin(agreement),"percolator_psm_score"].tolist()
    list_intersection_varcont=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(agreement),"percolator_psm_score"].tolist()
    list_pgonly=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),"percolator_psm_score"].tolist()
    plots.plot_ib_scores(list_ibonly,list_pgonly,list_intersection_varcont,list_intersection_varfree)
    #to explore: return scan ids and check the ids that were not identified with variant free method in a later function
    return(0)
    
def combidict_analysis(combidict,chromdict,stranddict,cpdt_pep,full_seqs,theoretical_saavs,mut_pep_probs,mut_cpdt_theoretical,mut_cpdt_counterparts,isOpenmut):
    # proteins_covered=Counter() #proteins detected
    mutated=set() #all proteins that were detected to have a variant by ionbot. how does compare to the proteins that actually do have variant?
    ref_only=set() #scan ids in the reference set
    ont_only=set() #scan ids in the ont set
    both=set() #scan ids that matched to both ref and ont proteins
    # hits_missed=0
    # hit_mut=0
    # hits_missed_mut=0
    chrom_dist=Counter() # 1 chromosome location per scan id
    strand_dist=Counter()
    mut_cpdt_observed=mut_cpdt_theoretical
    mut_counterparts_observed=mut_cpdt_counterparts
    protein_support=Counter()
    unamb_protsupport=Counter()
    for row in tqdm(combidict.iterrows()):
        scanid=row[1][0]
        mod=str(row[1][7])
        aamod=re.findall('[A-Z]->[A-Z]',mod)
        pep=row[1][3]
        prot_ids=row[1][9]
        if '||' in prot_ids:
            ids=prot_ids.split('||')
            for i in ids:
                protein_support[helper_functions.get_id(i)]+=1
        else: #unambiguous assignment!
            ids=[prot_ids]
            unamb_protsupport[helper_functions.get_id(prot_ids)]+=1
        # proteins_covered=detected_proteins(ids,proteins_covered) #what proteins from the proteome are covered and in what amounts
        cpdt_pep=helper_functions.add_to_observed(pep,ids,cpdt_pep,False) #what peptides are detected, how many, and what proteins they come from
        mut_counterparts_observed=helper_functions.add_to_observed(pep,ids,mut_counterparts_observed,isOpenmut) #note the occurences of "reference" versions of the SAAV peptides
        chrom_origin=helper_functions.find_chrom(ids,chromdict)
        strand_origin=helper_functions.find_strand(ids,stranddict)
        chrom_dist[chrom_origin]+=1 #which chromosome does the peptide belong to
        strand_dist[strand_origin]+=1 #which strand does the peptide belong to
        ref_only,ont_only,both=helper_functions.bin_hits_by_source(scanid,ids,ref_only,ont_only,both, isOpenmut)
        if isOpenmut and len(aamod)>0: #if ib detects a mutated peptide (for open variant search only)
            # hit_mut+=1
            mut_cpdt_observed=helper_functions.add_to_observed(pep,ids,mut_cpdt_observed,isOpenmut)
            for i in ids: 
                mutated.add(helper_functions.get_id(i))
        # elif isOpenmut and detect_mut_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)!='': ##very strange scenario here!!##
        #     print(scanid)
        elif not isOpenmut: #check if mutant peptide if not open mutation settings
            mut_cpdt_observed=helper_functions.add_to_observed(pep,ids,mut_cpdt_observed,isOpenmut)
    #create the figures
    # print("number of hits with detected variant = " +str(hit_mut)+ " matched to "+str(len(mutated))+ " proteins.")
    # print("number of mutant peptides not matched to predicted mutant peptides = " +str(hits_missed_mut))
    #create checkpoint- save the results from above so that whole analysis does not need to be repeated to re-create the graphs
    if isOpenmut:
        # with open('checkpoint.openmut.json','w'):
        #     json.dump()
        plots.plot_mut_abundance(mut_cpdt_observed,mut_counterparts_observed,cpdt_pep,full_seqs,"mutant_abundance_varfree.png")
        plots.plot_mut_vs_nonmut(mut_cpdt_observed,mut_counterparts_observed,theoretical_saavs,mut_pep_probs,"_varfree.png")
        plots.plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varfree.png","vertical_coverage_varfree.png")
        plots.plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varfree.png",isOpenmut)
        plots.plot_support(protein_support,unamb_protsupport,'protein_evidence_varfree.png')
    else:
        plots.plot_mut_abundance(mut_cpdt_observed,mut_counterparts_observed,cpdt_pep,full_seqs,"mutant_abundance_varcont.png")
        plots.plot_mut_vs_nonmut(mut_cpdt_observed,mut_counterparts_observed,theoretical_saavs,mut_pep_probs,"_varcont.png")
        plots.plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varcont.png","vertical_coverage_varcont.png")
        plots.plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varcont.png",isOpenmut)
        plots.plot_support(protein_support,unamb_protsupport,'protein_evidence_varcont.png')
    if isOpenmut:
        return(helper_functions.count_muts(mut_cpdt_observed),mutated,chrom_dist,strand_dist)
    return(helper_functions.count_muts(mut_cpdt_observed),chrom_dist,strand_dist)

