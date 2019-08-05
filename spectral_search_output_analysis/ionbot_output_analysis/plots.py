#!/usr/bin/env python3

import matplotlib, re, os,sys, itertools, argparse
from matplotlib_venn import venn2,venn2_unweighted,venn3, venn3_unweighted
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from collections import Counter
import calculations
import helper_functions

# Set the visualization settings (maybe need to be adjusted for saving figures to file)
matplotlib.rcParams['axes.titlesize'] = 'xx-large'
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

def plot_target_decoy(df, save_as, score_name='Ionbot psm score', plot_title='Search result'):
    """
    Plot for a given search engine output the target/decoy score distributions,
    the relation between the q_values and the PSM scores and a PP plot between
    the target and the decoy distribution.

    Positional arguments:
    df - Pandas DataFrame containing all rank 1 PSMs
    score_col - Name of the column containing the search engine scores
    decoy_col - Name of the column marking decoy PSMs as True and target PSMs as
    False
    qval_col - Name of the column containing each PSMs q-value

    Keyword arguments:
    score_name - score name used in axis labels
    plot_title - Plot title
    save_as - If not None, but string, save file to this filename
    """

    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    # Score distribution plot
    score_cutoff = df[(df['q_value'] <= 0.01) & (~df['DB'])].sort_values('q_value').iloc[-1]['ionbot_psm_score']
    plot_list = [list(x) for x in [df[df['DB']]['ionbot_psm_score'], df[~df['DB']]['ionbot_psm_score']]]
    axes[0].hist(plot_list, bins=30, label=['Decoy', 'Target'], color=['r', 'blue'], lw=1, rwidth=1)
    axes[0].vlines(x=score_cutoff, ymin=0, ymax=axes[0].get_ylim()[1], linestyles='dashed')
    axes[0].legend()
    axes[0].set_ylabel("Number of matches")
    axes[0].set_xlabel(score_name)
    #axes[0].set_xlim(0, 1)

    # Q value plot
    axes[1].plot(df.sort_values('ionbot_psm_score')['ionbot_psm_score'], df.sort_values('ionbot_psm_score')['q_value'])
    axes[1].vlines(x=score_cutoff, ymin=0, ymax=axes[1].get_ylim()[1], linestyles='dashed')
    axes[1].set_ylabel('q-value')
    axes[1].set_xlabel(score_name)
    #axes[1].set_xlim(0, 1)

    # PP plot
    ratio = df['DB'].value_counts()['D'] / df['DB'].value_counts()['T']
    Ft = ECDF(df[~df['DB']]['ionbot_psm_score'])
    Fd = ECDF(df[df['DB']]['ionbot_psm_score'])
    x = df[~df['DB']]['ionbot_psm_score']
    Fdp = Fd(x)
    Ftp = Ft(x)
    axes[2].scatter(Fdp, Ftp, s=4)
    axes[2].plot((0, 1), (0, ratio), color='r')
    axes[2].set_xlabel('Decoy percentile')
    axes[2].set_ylabel('Target percentile')

    plt.suptitle(plot_title)
    sns.despine()
    plt.savefig(save_as, facecolor='white', transparent=False)


def plot_qvalues_comparison(df_dict, q_value_col='q_value', decoy_col='DB', fdr_levels=None, title=''):
    """
    Plot number of identifications at all q-values for multiple datasets.

    df_dict: dict of `name -> dataframe`, with each dataframe containing
    q-values and booleans indicating whether the q-value belongs to
    a target or decoy PSM.
    q_value_col: Name of q-value column
    decoy_col: Name of decoy column
    fdr_levels: List of FDR float values to plot as vertical lines
    """
    d={'T':False,'D':True}
    for label, df_in in df_dict.items():
        df = df_in.reset_index(drop=True).sort_values(q_value_col, ascending=True).copy()
        df[decoy_col]=df[decoy_col].map(d)
        df['count'] = (~df[decoy_col]).cumsum()
        plt.plot(df[q_value_col], df['count'], label=label, alpha=0.5)

    for fdr in fdr_levels:
	    plt.plot([fdr]*2, np.linspace(0, np.max(df['count']), 2), linestyle='--', color='black', label='{} FDR'.format(fdr))
    plt.ylabel('Number of identified spectra')
    plt.xlabel('FDR (log scale)')
    #plt.xscale("log", nonposy='clip')
    plt.xscale("log")
    plt.xlim(0.00001, 1)
    plt.legend()
    plt.title(title)
    plt.savefig('qval_comparison.png')

def plot_support(prot_evidence,unamb_prot_evidence,figname):
    '''look into the support for proteins
    how many proteins have 1,2,3... peptides supporting their existence? unambiguously?
    1 counter that counts all peptides, another counter that counts only unique peptides'''
    plt.figure('support')
    recountpev=helper_functions.counter_translator(prot_evidence)
    recountunamb=helper_functions.counter_translator(unamb_prot_evidence)
    print("count #proteins with higher than 20 peptide evidence:"+str(recountpev['20+']))
    new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,'20+']
    pev= pd.DataFrame.from_dict(recountpev,orient='index')
    punamb= pd.DataFrame.from_dict(recountunamb,orient='index')
    combi=pd.concat([pev,punamb],axis=1,sort=True)
    combi=combi.reindex(new_index)
    combi.fillna(0,inplace=True)
    combi.columns=['Any peptide evidence','Unambiguous assignments']
    combi.plot(kind='bar',legend=False,title="Peptide support for proteins")
    plt.ylabel("# proteins")
    plt.ylim(0,20000)
    plt.xlabel("# peptides")
    plt.legend(loc='upper right')
    plt.savefig(figname)
    plt.clf()
    return('Plotted peptide support')

def plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_combi):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution")
    sns.distplot(ibdf_refonly['ionbot_psm_score'], hist=False, label='Human reference only',axlabel='Ionbot score')
    sns.distplot(ibdf_ontonly['ionbot_psm_score'], hist=False, label='Transcriptome translation only',axlabel='Ionbot score')
    sns.distplot(ibdf_combi['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation',axlabel='Ionbot score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptides in variant-free libraries')
    plt.savefig("qc_pearsonr_3source.png")
    plt.close()
    return("Scores plot made")

def plot_scores_pg(ibdf_combi,ibdf_combi_pg):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='T']['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation (no variants)',axlabel='Ionbot score')
    sns.distplot(ibdf_combi_pg[ibdf_combi_pg['DB']=='T']['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation (with variants)',axlabel='Ionbot score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr_pgvsopenmut.png")
    plt.close()
    return("Scores plot made")

def plot_scores_decoy(ibdf_combi,figname):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='T']['ionbot_psm_score'], label='Target',axlabel='Ionbot score')
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='D']['ionbot_psm_score'], label='Decoy',axlabel='Ionbot score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig(figname)
    plt.close()
    return("Scores plot made")

def plot_source_piechart(ref_only,ont_only,both,figname,isOpenmut):
    '''this function will plot the source piechart of sources of the hits and save it to a pdf'''
    plt.figure('source piechart')
    explode = (0.1, 0, 0)
    if isOpenmut:
        labels='Exclusively ONT transcriptome','Exclusively reference (gencode)', 'Both'
    else:
        labels='Novel','Annotated','Combination'
    plt.pie([len(ont_only),len(ref_only),len(both)],autopct='%1.1f%%', explode=explode,colors=['#de2d26','#3182bd','#756bb1'])
    plt.title('Peptide spectral hits by source',fontsize=35)
    plt.legend(labels,loc=8)
    plt.savefig(figname)
    plt.close()
    return("saved to sources_spectral_hits")

def plot_chromosomal_dist(distr_classic,distr_openmut):
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M','unknown']
    plt.figure('chromosomal distribution')
    chist=pd.DataFrame.from_dict(distr_classic,orient='index')#.sort_index()
    chist_openmut=pd.DataFrame.from_dict(distr_openmut,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_openmut],axis=1)
    combi=combi.reindex(new_index)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi.plot(kind='bar',legend=False,title="Chromosomal distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Chromosome")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('chromosomal_distribution.png')
    plt.close()
    return("plotted chromosomal distribution")

def plot_strand_dist(distr_classic,distr_openmut):
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('strand distribution')
    chist=pd.DataFrame.from_dict(distr_classic,orient='index')#.sort_index()
    chist_openmut=pd.DataFrame.from_dict(distr_openmut,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_openmut],axis=1)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi.plot(kind='bar',legend=False,title="Strand distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Strand")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('strand_distribution.png')
    plt.close()
    return("plotted strand distribution")

def plot_coverage_plots(cpdt_pep,fullseqs,fignamehorizontal,fignamevertical):
    '''this function will plot the graphs that correspond to the coverage of the proteome
    - vertical coverage
    - horizontal coverage
    - chromosome distribution
    '''
    high_cov_hor,cov_vert,perc_cov_dist=calculations.coverage_measure(cpdt_pep,fullseqs)
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    #horizontal coverage
    plt.figure('horizontal coverage')
    plt.hist(perc_cov_dist,bins=300)
    plt.xlim(1,100)
    plt.ylim(0,1000)
    plt.title("Horizontal coverage")
    plt.xlabel("% Coverage")
    plt.ylabel("# Proteins")
    plt.savefig(fignamehorizontal)
    plt.clf()
    #vertical coverage
    plt.figure('vertical coverage')
    plt.hist(cov_vert,bins=1000)
    plt.xlim(1,40)
    plt.ylim(0,4000)
    plt.title("Vertical coverage")
    plt.xlabel("Peptide count (protein size normalized)")
    plt.ylabel("Density")
    plt.savefig(fignamevertical)
    plt.close()
    return("Plotted coverage")

def plot_heatmaps(counter,prefix,suffix):
    '''plot the types of substitutions that occur'''
    outfile=prefix+suffix
    if type(counter)==Counter:
        ser=pd.Series(list(counter.values()),index=pd.MultiIndex.from_tuples(counter.keys()))
        df=ser.unstack()
    else:
        df=counter
    plt.figure("heatmap")
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    sns.heatmap(df)
    # plt.title("Substitutions")
    plt.ylabel("Original")
    plt.xlabel("Variant")
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    return(0)

def plot_mut_abundance(mutant_cpdtpep,counterpart_cpdtpep,cpdtpep,fullseqs,figname):
    '''plot protein abundance vs number of detected mutant peptides'''
    print("variant peptides:")
    mut_pep_abundance,nonmut_pep_abundance=calculations.calc_mut_abundances(mutant_cpdtpep,cpdtpep,fullseqs)
    print("non-variant counterparts:")
    cpt_pep_abundance,nmpa=calculations.calc_mut_abundances(counterpart_cpdtpep,cpdtpep,fullseqs)
    calculations.calculate_correlation(mut_pep_abundance,cpt_pep_abundance,nonmut_pep_abundance)
    #make plot
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('mutant peptides')
    # plt.scatter(*zip(*mut_pep_abundance),c='r',label='Variant peptide',alpha=1)
    # plt.scatter(*zip(*nonmut_pep_abundance),c='b',label='Normal peptide',alpha=0.25)
    x,y=zip(*mut_pep_abundance)
    sns.regplot(x=list(x),y=list(y),dropna=False,scatter=True,fit_reg=True, color='r',label='Variant peptide')
    x,y=zip(*nonmut_pep_abundance)
    sns.regplot(x=list(x),y=list(y),dropna=False,scatter=True,fit_reg=True,color='b',label='Normal peptide')
    plt.xlabel('Protein abundance (NSAF normalized)')
    plt.xlim(0,0.0002)
    plt.ylabel('Number mutant peptides detected')
    plt.ylim(-10,1000)
    plt.title('Peptide abundance vs total protein abundance')
    plt.legend(loc='upper right')
    plt.savefig(figname)
    plt.close()
    return('done')

def plot_mut_vs_prob(counts,figname):
    '''plot variant observed count vs CPDT probability for that peptide'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('measure probability vs observed peptides')
    sns.jointplot(*zip(*counts), kind="reg", stat_func=calculations.r2)
    plt.xlabel('Variant peptide theoretical detectability')
    plt.ylabel('Variant peptide observed count')
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    return(0)

def plot_mut_vs_nonmut(mutant_cpdtpep,counterpart_cpdtpep,theoretical_counts,var_probs,suffix):
    counts,probs,observed=calculations.calc_pep_counts(mutant_cpdtpep,counterpart_cpdtpep,var_probs)
    plot_heatmaps(observed,'heatmap_observed_subs',suffix)
    comparison=helper_functions.get_normalized_matrix(theoretical_counts,observed)
    plot_heatmaps(comparison,'heatmap_difference',suffix)
    plot_mut_vs_prob(probs,'probability_detection_vs_count'+suffix)
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('measure direct counterparts')
    # sns.regplot(*zip(*counts),scatter=True,fit_reg=True,color='b',alpha=1)
    sns.jointplot(*zip(*counts), kind="reg", stat_func=calculations.r2)
    plt.xlabel('Variant peptide count')
    plt.ylabel('Reference counterpart count')
    plt.ylim(-10,700)
    # plt.title('Variant vs. non-variant peptide abundance')
    plt.tight_layout()
    plt.savefig("variant_vs_nonvariant"+suffix)
    plt.close()
    return(0)

def plot_ib_scores_directcomp(varfree_scores,varcont_scores):
    '''for the variant peptides that were found in the variant containing set but not in the variant free set,
    what is the ionbot score distribution from each respective results list'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("ionbot scores discrepant hits")
    #inner join the 2
    combi=pd.merge(varfree_scores,varcont_scores,on="scan_id",suffixes=("_varfree","_varcont"))
    combi.groupby("matched_peptide").mean() #make sure don't have groups of dots per unique peptide
    combi["pep_length"]=combi["matched_peptide"].str.len() #record length of matched peptide (by var-free)
    combi.plot.scatter(x="ionbot_psm_score_varfree",y="ionbot_psm_score_varcont",c="pep_length",colormap='viridis')
    # sns.distplot(varfree_scores, hist=False, label='Variant-free',axlabel='Ionbot score')
    # sns.distplot(varcont_scores, hist=False, label='Variant-containing',axlabel='Ionbot score')
    plt.ylabel("Scores variant-containing")
    plt.xlabel("Scores variant-free")
    plt.legend()
    plt.title('Ionbot scores for discrepant peptides found only in the variant-containing search')
    plt.savefig("discrepant_peptide_direct_comparison.png")
    plt.close()
    return("Scores plot made")

def plot_ib_scores(ibonly,pgonly,intersectionpg,intersectionom):
    '''for the variant peptides that were found in the variant containing set but not in the variant free set,
    what is the ionbot score distribution from each respective results list'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("ionbot scores discrepant hits")
    # sns.distplot(ibonly, hist=False, label='Variant-free only',axlabel='Ionbot score')
    # sns.distplot(pgonly, hist=False, label='Variant-containing only',axlabel='Ionbot score')
    # sns.distplot(intersectionpg, hist=False, label='Intersection variant-containing',axlabel='Ionbot score')
    # if len(intersectionom)>0:
    #     sns.distplot(intersectionom, hist=False, label='Intersection variant-free')
    plt.boxplot([ibonly,pgonly,intersectionpg,intersectionom])
    plt.xticks([1,2,3,4],['Variant-free only','Variant-containing only','Intersection variant-containing','Intersection variant-free'])
    # plt.legend()
    plt.title('Ionbot scores for detected variant peptides')
    plt.savefig("discrepant_peptide_scores.png")
    plt.close()
    return("Scores plot made")

def plot_unexpected_mods(mods_om,mods_pg):
    mod_ct_om=helper_functions.categorize_mods(mods_om)
    mod_ct_pg=helper_functions.categorize_mods(mods_pg)
    plt.figure('discrepant peptide lengths')
    chist_pg=pd.DataFrame.from_dict(dict(mod_ct_pg.most_common(10)),orient='index')
    chist_om=pd.DataFrame.from_dict(dict(mod_ct_om.most_common(10)),orient='index')
    combi=pd.concat([chist_pg,chist_om],axis=1,sort=False)
    combi.fillna(0)
    combi.columns=['Variant-containing','Variant-free']
    combi.plot(kind='bar',title="Unexpected modifications found instead of SAAVs from variant peptides (variant-free search)")
    plt.ylabel("Count peptides")
    plt.xlabel("PTMs")
    plt.legend()
    plt.tight_layout()
    plt.savefig('discrepant_peptide_mods.png')
    plt.close()
    return(0)

def plot_peplengths(peptide_counter_pg,peptide_counter_om):
    lenct_pg=helper_functions.gather_counts(peptide_counter_pg)
    lenct_om=helper_functions.gather_counts(peptide_counter_om)
    plt.figure('discrepant peptide lengths')
    # new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M','unknown']
    chist_pg=pd.DataFrame.from_dict(lenct_pg,orient='index',columns=["Combi variant-containing only"]).sort_index()
    if len(lenct_om)>0:
        chist_om=pd.DataFrame.from_dict(lenct_om,orient='index',columns=["Combi variant-free only"]).sort_index()
        combi=pd.concat([chist_pg,chist_om],axis=1)
        combi.fillna(0)
    else:
        combi=chist_pg
        combi["Combi variant-free only"]=0
    combi.plot(kind='bar',legend=False,title="Length of discrepant peptides (found by one method and not the other)")
    # combi=combi.reindex(new_index)
    plt.ylabel("Density")
    plt.xlabel("Length peptide")
    plt.legend(loc='upper left')
    plt.savefig('discrepant_peptide_length.png')
    plt.close()
    return(0)

def plot_final_venns(allmuts_classic,allmuts_openmut,mut_cpdt_theoretical,mutprotset):
    #create diagrams
    plt.figure('venn mutant psms')
    vda=venn2_unweighted([allmuts_classic,allmuts_openmut],('Combi variant-containing','Combi variant-free')) #venn for the overlap in detected peptides
    plt.title("Observed variant PSMs",fontsize=26)
    for text in vda.set_labels:
        text.set_fontsize(26)
    for text in vda.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_peps.png')
    plt.clf()
    plt.figure('venn mutant peptides')
    vdb=venn2_unweighted([set(allmuts_classic),set(allmuts_openmut)],("Combi variant-containing","Combi variant-free")) #venn for the overlap in detected proteins
    plt.title("Unique observed variant peptides",fontsize=26)
    for text in vdb.set_labels:
        text.set_fontsize(26)
    for text in vdb.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_prots.png')
    plt.clf()
    plt.figure('venn proteins all')
    vdb=venn3([set(allmuts_classic),set(mut_cpdt_theoretical.keys()),mutprotset],("Combi variant-containing","All theoretical","Combi variant-free")) #venn for the overlap in detected proteins
    plt.title("Unique proteins associated with variant peptides",fontsize=26)
    # for text in vdb.set_labels:
    #     text.set_fontsize(26)
    # for text in vdb.subset_labels:
    #     text.set_fontsize(20)
    plt.savefig('overlap_all_detected_mut_prots.png')
    plt.close()
    return('plotted final venns')

