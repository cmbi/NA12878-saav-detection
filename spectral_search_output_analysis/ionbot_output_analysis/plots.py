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
import logging
log = logging.getLogger(__name__)

# Set the visualization settings (maybe need to be adjusted for saving figures to file)
# matplotlib.rcParams['axes.titlesize'] = 'xx-large'
# matplotlib.rcParams['axes.labelsize'] = 'x-large'
# matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
matplotlib.rcParams.update({'font.size': 22})

def plot_target_decoy(df, save_as, score_name='Percolator psm score', plot_title='Search result'):
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

    score_cutoff = df[(df['q_value'] <= 0.01) & (~df['DB'])].sort_values('q_value').iloc[-1]['percolator_psm_score']
    plot_list = [list(x) for x in [df[df['DB']]['percolator_psm_score'], df[~df['DB']]['percolator_psm_score']]]
    axes[0].hist(plot_list, bins=30, label=['Decoy', 'Target'], color=['r', 'blue'], lw=1, rwidth=1)
    axes[0].vlines(x=score_cutoff, ymin=0, ymax=axes[0].get_ylim()[1], linestyles='dashed')
    axes[0].legend()
    axes[0].set_ylabel("Number of matches")
    axes[0].set_xlabel(score_name)
    #axes[0].set_xlim(0, 1)

    # Q value plot
    axes[1].plot(df.sort_values('percolator_psm_score')['percolator_psm_score'], df.sort_values('percolator_psm_score')['q_value'])
    axes[1].vlines(x=score_cutoff, ymin=0, ymax=axes[1].get_ylim()[1], linestyles='dashed')
    axes[1].set_ylabel('q-value')
    axes[1].set_xlabel(score_name)
    #axes[1].set_xlim(0, 1)

    # PP plot
    ratio = df['DB'].value_counts()[True] / df['DB'].value_counts()[False]
    Ft = ECDF(df[~df['DB']]['percolator_psm_score'])
    Fd = ECDF(df[df['DB']]['percolator_psm_score'])
    x = df[~df['DB']]['percolator_psm_score']
    Fdp = Fd(x)
    Ftp = Ft(x)
    axes[2].scatter(Fdp, Ftp, s=4)
    axes[2].plot((0, 1), (0, ratio), color='r')
    axes[2].set_xlabel('Decoy percentile')
    axes[2].set_ylabel('Target percentile')

    plt.suptitle(plot_title)
    plt.tight_layout()
    sns.despine()
    plt.savefig(save_as, facecolor='white', transparent=False)
    plt.close()


def plot_qvalues_comparison(df_dict, q_value_col='q_value', decoy_col='DB', fdr_levels=None):
    """
    Plot number of identifications at all q-values for multiple datasets.

    df_dict: dict of `name -> dataframe`, with each dataframe containing
    q-values and booleans indicating whether the q-value belongs to
    a target or decoy PSM.
    q_value_col: Name of q-value column
    decoy_col: Name of decoy column
    fdr_levels: List of FDR float values to plot as vertical lines
    """
    matplotlib.rcParams.update({'font.size': 22})

    for label, df_in in df_dict.items():
        df = df_in.reset_index(drop=True).sort_values(q_value_col, ascending=True).copy()
        df[decoy_col]=df[decoy_col]=='D'
        df['count'] = (~df[decoy_col]).cumsum()
        plt.plot(df[q_value_col], df['count'], label=label, alpha=0.5)

    for fdr in fdr_levels:
	    plt.plot([fdr]*2, np.linspace(0, np.max(df['count']), 2), linestyle='--', color='black', label='{} FDR'.format(fdr))
    plt.ylabel('Number of identified spectra')
    plt.xlabel('FDR (log scale)')
    #plt.xscale("log", nonposy='clip')
    plt.xscale("log")
    plt.xlim(0.0001, 1)
    plt.legend()
    plt.title('q-value comparison search dictionaries')
    plt.tight_layout()
    plt.savefig('qval_comparison.png')
    
    plt.close()

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
    
    plt.close()

def plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_vf):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution")
    sns.distplot(ibdf_refonly['percolator_psm_score'], hist=False, label='Human reference only',axlabel='Percolator score')
    sns.distplot(ibdf_ontonly['percolator_psm_score'], hist=False, label='Transcriptome translation only',axlabel='Percolator score')
    sns.distplot(ibdf_vf['percolator_psm_score'], hist=False, label='Reference + transcriptome translation',axlabel='Percolator score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptides in variant-free libraries')
    plt.savefig("qc_pearsonr_3source.png")
    
    plt.close()

def plot_scores_combi(ibdf_vf,ibdf_vc):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_vf[ibdf_vf['DB']=='T']['percolator_psm_score'], hist=False, label='Reference + transcriptome translation (no variants)',axlabel='Percolator score')
    sns.distplot(ibdf_vc[ibdf_vc['DB']=='T']['percolator_psm_score'], hist=False, label='Reference + transcriptome translation (with variants)',axlabel='Percolator score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr_vcvsopenmut.png")
    plt.close()

def plot_scores_decoy(ibdf_vf,figname):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_vf[ibdf_vf['DB']=='T']['percolator_psm_score'], label='Target',axlabel='Percolator score')
    sns.distplot(ibdf_vf[ibdf_vf['DB']=='D']['percolator_psm_score'], label='Decoy',axlabel='Percolator score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig(figname)
    
    plt.close()

def plot_source_piechart(source_counter,figname):
    '''this function will plot the source piechart of sources of the hits and save it to a pdf'''
    plt.figure('source piechart')
    explode = (0.1, 0, 0)
    labels='Exclusively ONT transcriptome','Exclusively reference (gencode)', 'Both'
    plt.pie([source_counter['ont'],source_counter['ref'],source_counter['both']],autopct='%1.1f%%', explode=explode,colors=['#de2d26','#3182bd','#756bb1'])
    # plt.title('Peptide spectral hits by source',fontsize=35)
    plt.legend(labels,loc=8)
    plt.savefig(figname)
    
    plt.close()

def plot_chromosomal_dist(distr_vc,distr_vf):
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    # new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M','unknown']
    plt.figure('chromosomal distribution')
    chist=pd.DataFrame.from_dict(distr_vc,orient='index')#.sort_index()
    chist_vf=pd.DataFrame.from_dict(distr_vf,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_vf],axis=1,sort=True)
    # combi=combi.reindex(new_index)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi.plot(kind='bar',legend=False,title="Chromosomal distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Chromosome")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('chromosomal_distribution.png')
    
    plt.close()

def plot_strand_dist(distr_vc,distr_vf):
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('strand distribution')
    chist=pd.DataFrame.from_dict(distr_vc,orient='index')#.sort_index()
    chist_vf=pd.DataFrame.from_dict(distr_vf,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_vf],axis=1)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi.plot(kind='bar',legend=False,title="Strand distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Strand")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('strand_distribution.png')
    
    plt.close()

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

def plot_heatmaps(counter,outfile):
    '''plot the types of substitutions that occur'''
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
    sns.regplot(x=list(x),y=list(y),dropna=False,scatter=True,fit_reg=True, color='b')
    # x,y=zip(*nonmut_pep_abundance)
    # sns.regplot(x=list(x),y=list(y),dropna=False,scatter=True,fit_reg=True,color='b',label='Normal peptide')
    plt.xlabel('Non-variant peptide abundance (NSAF normalized)')
    plt.xlim(0,0.0002)
    plt.ylabel('Count variant peptides')
    plt.ylim(-10,1000)
    plt.title('Variant peptide abundance vs non-variant peptide abundance')
    # plt.legend(loc='upper right')
    plt.savefig(figname)
    plt.close()
    

def plot_mut_vs_prob(counts,figname):
    '''plot variant observed count vs CPDT probability for that peptide'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('measure probability vs observed peptides')
    sns.jointplot(*zip(*counts))
    plt.xlabel('Variant peptide theoretical detectability')
    plt.ylabel('Variant peptide observed count')
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def plot_mut_vs_nonmut(counts,figname):
    # sns.set(style="white", color_codes=True)
    plt.figure('measure direct counterparts')
    # sns.regplot(*zip(*counts),scatter=True,fit_reg=True,color='b',alpha=1)
    h=sns.jointplot(*zip(*counts), kind='scatter',stat_func=calculations.r2)
    h.set_axis_labels('Variant peptide count', 'Reference counterpart count', fontsize=16)
    # ax = h.ax_joint
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # h.ax_marg_x.set_xlim(0, 250)
    h.ax_marg_x.set_xscale('log')
    h.ax_marg_y.set_yscale('log')
    left, right = plt.xlim()
    x = np.linspace(left,right)
    h.ax_joint.plot(x,x,':k')
    # plt.title('Variant vs. non-variant peptide abundance')
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def plot_ib_scores_directcomp(combi,retentiontime):
    '''for the variant peptides that were found in one set but not another,
    what is the Percolator score distribution from each respective results list
    color by retention time prediction instead of length
    '''
    rt=pd.read_csv(retentiontime)
    rt.columns=['matched_peptide','predicted_retention_time']
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Percolator scores discrepant hits")
    # combi.groupby("matched_peptide").mean()
    # combi.groupby("matched_peptide_vf").mean() #make sure don't have groups of dots per unique peptide
    # combi["pep_length"]=combi["matched_peptide"].str.len() #record length of matched peptide (by var-free)
    combi=pd.merge(combi,rt,on='matched_peptide')
    # combi.plot.scatter(x="percolator_psm_score_varfree",y="percolator_psm_score_varcont",c="pep_length",colormap='viridis')
    combi.plot.scatter(x="percolator_psm_score_vf",y="percolator_psm_score_vc",c="predicted_retention_time",colormap='viridis')
    left, right = plt.xlim()
    x = np.linspace(left,right)
    plt.plot(x, x)
    # sns.distplot(varfree_scores, hist=False, label='Variant-free',axlabel='Percolator score')
    # sns.distplot(varcont_scores, hist=False, label='Variant-containing',axlabel='Percolator score')
    plt.ylabel("Scores variant-containing")
    plt.xlabel("Scores variant-free")
    plt.legend()
    plt.title('Percolator scores for discrepant peptides found only in the variant-containing search')
    plt.savefig("discrepant_peptide_direct_comparison.png")
    
    plt.close()

def plot_ib_scores(ibonly,pgonly,intersectionpg,intersectionom,nonmutvc,nonmutvf):
    '''for the variant peptides that were found in the variant containing set but not in the variant free set,
    what is the Percolator score distribution from each respective results list'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Percolator scores discrepant hits")
    plt.boxplot([ibonly,pgonly,intersectionpg,intersectionom,nonmutvc,nonmutvf])
    plt.xticks([1,2,3,4,5,6],['VF only','VC only','Intersection VC','Intersection VF','Non-variant VC','Non-variant VF'])
    # plt.legend()
    plt.title('Percolator scores for detected variant peptides')
    plt.savefig("discrepant_peptide_scores.png")
    plt.clf()
    plt.close()

def plot_unexpected_mods(mods_vf,mods_vc,variant=False):
    mod_ct_vf=Counter(dict(mods_vf.value_counts()))
    plt.figure('discrepant peptide lengths')
    chist_vf=pd.DataFrame.from_dict(dict(mod_ct_vf.most_common(10)),orient='index')
    if variant:
        figname='variant_unexpected_mods.png'
        title="Unexpected modifications found instead of SAAVs from variant peptides (variant-free search)"
        chist_vf.plot(kind='bar',title=title,legend=False)
    else:
        figname='nonvariant_unexpected_mods.png'
        title='All unexpected modifications found in non-variant peptides'
        mod_ct_vc=Counter(dict(mods_vc.value_counts()))
        chist_vc=pd.DataFrame.from_dict(dict(mod_ct_vc.most_common(10)),orient='index')
        combi=pd.concat([chist_vc,chist_vf],axis=1,sort=False)
        combi.fillna(0)
        combi.columns=['Nonvariant VC', 'Nonvariant VF']
        combi.plot(kind='bar',title=title)
    plt.ylabel("Count peptides")
    plt.xlabel("PTMs")
    plt.tight_layout()
    plt.savefig(figname)
    plt.clf()
    plt.close()

def plot_peplengths(lenct_vc,lenct_vf,variant=False):
    if variant:
        labels=['Variant VC','Variant VF']
        figname='variant_peptide_length.png'
    else:
        labels=['Nonvariant VC','Nonvariant VF']
        figname='nonvariant_peptide_length.png'
    plt.figure('discrepant peptide lengths')
    # new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M','unknown']
    chist_vc=pd.DataFrame.from_dict(lenct_vc,orient='index').sort_index()
    if len(lenct_vf)>0:
        chist_vf=pd.DataFrame.from_dict(lenct_vf,orient='index').sort_index()
        combi=pd.concat([chist_vc,chist_vf],axis=1)
        combi.fillna(0)
    else:
        combi=chist_vc
        labels=['Variant VC']
    combi.columns=labels
    combi.plot(kind='bar',legend=False,title="Length of variant and normal peptides from combination search dictionaries")
    # combi=combi.reindex(new_index)
    plt.ylabel("Density")
    plt.xlabel("Length peptide")
    plt.legend(loc='upper left')
    plt.savefig(figname)
    plt.clf()
    plt.close()

def plot_final_venns(df_vc,df_vf):
    #get appropriate counters
    vc=Counter(dict(df_vc['peptide'].value_counts()))
    vf=Counter(dict(df_vf['peptide'].value_counts()))
    #create diagrams
    plt.figure('venn variant psms')
    vda=venn2([vc,vf],('Combi variant-containing','Combi variant-free')) #venn for the overlap in detected peptides
    plt.title("Observed variant PSMs",fontsize=26)
    for text in vda.set_labels:
        text.set_fontsize(26)
    for text in vda.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_psms.png')
    plt.clf()
    plt.figure('venn variant peptides')
    vdb=venn2([set(vc),set(vf)],("Combi variant-containing","Combi variant-free")) #venn for the overlap in detected proteins
    plt.title("Unique observed variant peptides",fontsize=26)
    for text in vdb.set_labels:
        text.set_fontsize(26)
    for text in vdb.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_peps.png')
    plt.clf()
    # plt.figure('venn proteins all')
    # vdb=venn3([set(vc),set(mut_cpdt_theoretical.keys()),mutprotset],("Combi variant-containing","All theoretical","Combi variant-free")) #venn for the overlap in detected proteins
    # plt.title("Unique proteins associated with variant peptides",fontsize=26)
    # # for text in vdb.set_labels:
    # #     text.set_fontsize(26)
    # # for text in vdb.subset_labels:
    # #     text.set_fontsize(20)
    # plt.savefig('overlap_all_detected_mut_prots.png')
    plt.close()

