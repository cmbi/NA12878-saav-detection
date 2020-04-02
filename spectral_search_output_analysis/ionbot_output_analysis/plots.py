#!/usr/bin/env python3

import matplotlib, re, os,sys, itertools, argparse
from matplotlib_venn import venn2,venn2_unweighted,venn3, venn3_unweighted
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF
from collections import Counter
import calculations
import helper_functions
import logging
# log = logging.getLogger(__name__)

# Set the visualization settings (maybe need to be adjusted for saving figures to file)
# matplotlib.rcParams['axes.titlesize'] = 'xx-large'
# matplotlib.rcParams['axes.labelsize'] = 'x-large'
# matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
# matplotlib.rcParams.update({'font.size': 22})
sns.set_context('paper')

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

    score_cutoff = df[(df['q_value'] <= 0.01) & (~df['DB'])].sort_values('q_value').iloc[-1]['percolator_psm_score_best']
    plot_list = [list(x) for x in [df[df['DB']]['percolator_psm_score_best'], df[~df['DB']]['percolator_psm_score_best']]]
    axes[0].hist(plot_list, bins=30, label=['Decoy', 'Target'], color=['r', 'blue'], lw=1, rwidth=1)
    axes[0].vlines(x=score_cutoff, ymin=0, ymax=axes[0].get_ylim()[1], linestyles='dashed')
    axes[0].legend()
    axes[0].set_ylabel("Number of matches")
    axes[0].set_xlabel(score_name)
    #axes[0].set_xlim(0, 1)

    # Q value plot
    axes[1].plot(df.sort_values('percolator_psm_score_best')['percolator_psm_score_best'], df.sort_values('percolator_psm_score_best')['q_value'])
    axes[1].vlines(x=score_cutoff, ymin=0, ymax=axes[1].get_ylim()[1], linestyles='dashed')
    axes[1].set_ylabel('q-value')
    axes[1].set_xlabel(score_name)
    #axes[1].set_xlim(0, 1)

    # PP plot
    ratio = df['DB'].value_counts()[True] / df['DB'].value_counts()[False]
    Ft = ECDF(df[~df['DB']]['percolator_psm_score_best'])
    Fd = ECDF(df[df['DB']]['percolator_psm_score_best'])
    x = df[~df['DB']]['percolator_psm_score_best']
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


def plot_qvalues_comparison(df_dict, q_value_col='q_value', decoy_col='DB', fdr_levels=None,plotname='qval_comparison.png'):
    """
    Plot number of identifications at all q-values for multiple datasets.

    df_dict: dict of `name -> dataframe`, with each dataframe containing
    q-values and booleans indicating whether the q-value belongs to
    a target or decoy PSM.
    q_value_col: Name of q-value column
    decoy_col: Name of decoy column
    fdr_levels: List of FDR float values to plot as vertical lines
    """
    #matplotlib.rcParams.update({'font.size': 22})
    for label, df_in in df_dict.items():
        df = df_in.reset_index(drop=True).sort_values(q_value_col, ascending=True).copy()
        df['count'] = (~df[decoy_col]).cumsum()
        plt.plot(df[q_value_col], df['count'], label=label, alpha=0.5)
    for fdr in fdr_levels:
	    plt.plot([fdr]*2, np.linspace(0, np.max(df['count']), 2), linestyle='--', color='black', label='{} FDR'.format(fdr))
    plt.ylabel('Number of identified spectra')
    plt.xlabel('FDR (log scale)')
    #plt.xscale("log", nonposy='clip')
    plt.xscale("log")
    plt.xlim(0.0001, 1)
    plt.legend(loc=2,fontsize='xx-small')
    plt.title('q-value comparison search dictionaries')
    plt.tight_layout()
    plt.savefig(plotname)
    
    plt.close()

def plot_support(prot_evidence,unamb_prot_evidence,figname):
    '''look into the support for proteins
    how many proteins have 1,2,3... peptides supporting their existence? unambiguously?
    1 counter that counts all peptides, another counter that counts only unique peptides'''
    plt.figure('support')
    #sns.set(font_scale=2)
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
    #plt.ylim(0,20000)
    plt.xlabel("# peptides")
    plt.legend(loc='upper right')
    plt.savefig(figname)
    
    plt.close()

def plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_vf):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    #sns.set(rc={'figure.figsize':(11.7,8.27)})
    #sns.set_style(style='white')
    sns.set_context('paper')
    plt.figure("Pearson R distribution")
    sns.distplot(ibdf_refonly['percolator_psm_score_best'], hist=False, label='Human reference only',axlabel='Percolator score')
    sns.distplot(ibdf_ontonly['percolator_psm_score_best'], hist=False, label='Transcriptome translation only',axlabel='Percolator score')
    sns.distplot(ibdf_vf['percolator_psm_score_best'], hist=False, label='Combi variant-free',axlabel='Percolator score')
    plt.legend(loc=1)
    plt.title('Correlation between theoretical and observed PSMs in variant-free libraries')
    plt.savefig("qc_pearsonr_3source.png")
    
    plt.close()

def quickplot(dfdict,filetype):
    '''plot the score distribution of one of the csv files'''
    sns.set_context('paper')
    fig = plt.figure()
    ax = plt.subplot(111)
    fig.set_size_inches(18.5,10.5)
    for title,df in dfdict.items():
        sns.distplot(df['percolator_psm_score_best'], hist=False, label=title,axlabel='Percolator score')
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize='xx-small')
    #plt.legend(loc=1)
    plt.savefig(f'{filetype}_dist_files_qc.png')
    plt.close()

def plot_scores_combi(ibdf_vf,ibdf_vc):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    #sns.set(rc={'figure.figsize':(11.7,8.27)},font_scale=2)
    #sns.set_style(style='white')
    sns.set_context('paper')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_vf[ibdf_vf['DB']==False]['percolator_psm_score_best'], hist=False, label='Combi variant-free',axlabel='Percolator score')
    sns.distplot(ibdf_vc[ibdf_vc['DB']==False]['percolator_psm_score_best'], hist=False, label='Combi variant-containing',axlabel='Percolator score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_score_vcvsopenmut.png")
    plt.close()

def plot_scores_decoy(ibdf_vf,figname):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    #sns.set(rc={'figure.figsize':(11.7,8.27)},font_scale=2)
    #sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_vf[ibdf_vf['DB']=='T']['percolator_psm_score_best'], label='Target',axlabel='Percolator score')
    sns.distplot(ibdf_vf[ibdf_vf['DB']=='D']['percolator_psm_score_best'], label='Decoy',axlabel='Percolator score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig(figname)
    
    plt.close()

def plot_source_piechart(source_counter,figname):
    '''this function will plot the source piechart of sources of the hits and save it to a pdf'''
    plt.figure('source piechart')
    #sns.set(rc={'figure.figsize':(11.7,8.27)},font_scale=2)
    #sns.set_style(style='white')
    explode = (0.1, 0, 0)
    labels='Exclusively ONT transcriptome','Exclusively reference (gencode)', 'Both'
    plt.pie([source_counter['ont'],source_counter['ref'],source_counter['both']],autopct='%1.1f%%', explode=explode,colors=['#de2d26','#3182bd','#756bb1'])
    # plt.title('Peptide spectral hits by source',fontsize=35)
    plt.legend(labels,loc=8)
    plt.savefig(figname)
    
    plt.close()

def plot_chromosomal_dist(distr_vc,distr_vf):
    #sns.set(rc={'figure.figsize':(11.7,8.27)},font_scale=2)
    #sns.set_style(style='white')
    new_index= ['chr1', 'chr2', 'chr3', 'chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    plt.figure('chromosomal distribution')
    chist=pd.DataFrame.from_dict(distr_vc,orient='index')#.sort_index()
    chist_vf=pd.DataFrame.from_dict(distr_vf,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_vf],axis=1,sort=True)
    combi=combi.reindex(new_index)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi['Combi variant-containing']=combi['Combi variant-containing']*100
    combi['Combi variant-free']=combi['Combi variant-free']*100
    combi.plot(kind='bar',legend=False,title="Chromosomal distribution of peptide hits")
    plt.ylabel("%  Identified peptides")
    plt.xlabel("Chromosome")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('chromosomal_distribution.png')
    
    plt.close()

def plot_strand_dist(distr_vc,distr_vf):
    # sns.set(rc={'figure.figsize':(11.7,8.27)},font_scale=2)
    # sns.set_style(style='white')
    plt.figure('strand distribution')
    chist=pd.DataFrame.from_dict(distr_vc,orient='index')#.sort_index()
    chist_vf=pd.DataFrame.from_dict(distr_vf,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_vf],axis=1,sort=True)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi['Combi variant-containing']=combi['Combi variant-containing']*100
    combi['Combi variant-free']=combi['Combi variant-free']*100
    combi.plot(kind='bar',legend=False,title="Strand distribution of peptide hits")
    plt.ylabel("%  Identified peptides")
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
    # sns.set(rc={'figure.figsize':(11.7,8.27)})
    # sns.set_style(style='white')
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

def plot_heatmaps(counter,outfile,bl=False):
    '''plot the types of substitutions that occur'''
    if type(counter)==dict:
        counter=pd.DataFrame.from_dict(counter,orient='index')
        counter.columns=['substitution']
        #matrix=ser.unstack()
    elif type(counter)==pd.Series:
        counter=counter.apply(lambda x: tuple(x.split(','))).value_counts()
    df=helper_functions.initiate_counter(bl).merge(pd.DataFrame(counter),left_on='sub',right_index=True,how='left').fillna(0).set_index('sub')
    df.index=pd.MultiIndex.from_tuples(df.index,names=('original','new'))
    matrix=df.unstack()#.to_numpy()
    plt.figure("heatmap")
    sns.set(rc={'figure.figsize':(11.7,11)},font_scale=2)
    sns.set_style(style='white')
    sns.heatmap(matrix['substitution'])
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
    sns.set_context('paper')
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
    

# def plot_mut_vs_prob(counts,figname):
#     '''plot variant observed count vs CPDT probability for that peptide'''
#     sns.set(rc={'figure.figsize':(11.7,8.27)})
#     sns.set_style(style='white')
#     plt.figure('measure probability vs observed peptides')
#     sns.jointplot(*zip(*counts))
#     plt.xlabel('Variant peptide theoretical detectability')
#     plt.ylabel('Variant peptide observed count')
#     plt.tight_layout()
#     plt.savefig(figname)
#     plt.close()

def plot_mut_vs_nonmut(counts,figname):
    # sns.set(style="white", color_codes=True)
    plt.figure('measure direct counterparts')
    # sns.regplot(*zip(*counts),scatter=True,fit_reg=True,color='b',alpha=1)
    het=counts[counts['var_type']=='heterozygous']
    hom=counts[counts['var_type']=='homozygous']
    h=sns.jointplot(het['count_var'],het['count_refctp'], kind='scatter',alpha=0.55,color='r',stat_func=calculations.r2)
    h.x=hom['count_var']
    h.y=hom['count_refctp']
    h.plot_joint(plt.scatter,alpha=0.55,color='b')
    h.set_axis_labels('Variant peptide count', 'Reference counterpart count', fontsize=16)
    # h.ax_marg_x.set_xlim(0, 250)
    #h.ax_marg_x.set_xscale('log')
    #h.ax_marg_y.set_yscale('log')
    left, right = plt.xlim()
    x = np.linspace(left,right)
    h.ax_joint.plot(x,x,':k')
    # plt.title('Variant vs. non-variant peptide abundance')
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def plot_ib_scores_directcomp(combi):
    '''for the variant peptides that were found in one set but not another,
    what is the Percolator score distribution from each respective results list
    color by retention time prediction instead of length
    '''
    sns.set_context('paper')
    plt.figure("Percolator scores discrepant hits")
    # combi.groupby("matched_peptide").mean()
    # combi.groupby("matched_peptide_vf").mean() #make sure don't have groups of dots per unique peptide
    # combi["pep_length"]=combi["matched_peptide"].str.len() #record length of matched peptide (by var-free)
    # combi=pd.merge(combi,rt_pred,on='matched_peptide')
    combi['delta_retention_time']=(combi['rt']-combi['predicted_tr']).abs()
    # combi.plot.scatter(x="percolator_psm_score_best_varfree",y="percolator_psm_score_best_varcont",c="pep_length",colormap='viridis')
    combi.plot.scatter(x="percolator_psm_score_best_vf",y="percolator_psm_score_best_vc",c="delta_retention_time",colormap='viridis')
    left, right = plt.xlim()
    x = np.linspace(left,right)
    plt.plot(x, x)
    # sns.distplot(varfree_scores, hist=False, label='Variant-free',axlabel='Percolator score')
    # sns.distplot(varcont_scores, hist=False, label='Variant-containing',axlabel='Percolator score')
    plt.ylabel("Scores variant-containing")
    plt.xlabel("Scores variant-free")
    plt.legend()
    plt.title('VF percolator scores for variant peptides found only in the VC search')
    plt.savefig("discrepant_peptide_direct_comparison.png")
    
    plt.close()

def plot_ib_scores(ibonly,pgonly,intersectionpg,intersectionom,nonmutvc,nonmutvf):
    '''for the variant peptides that were found in the variant containing set but not in the variant free set,
    what is the Percolator score distribution from each respective results list'''
    sns.set_context('paper')
    plt.figure("Percolator scores discrepant hits")
    plt.boxplot([ibonly,pgonly,intersectionpg,intersectionom,nonmutvc,nonmutvf])
    plt.xticks([1,2,3,4,5,6],['VF only','VC only','Intersection VC','Intersection VF','Non-variant VC','Non-variant VF'])
    # plt.legend()
    plt.title('Percolator scores for detected variant peptides')
    plt.savefig("discrepant_peptide_scores.png")
    plt.clf()
    plt.close()

def plot_unexpected_mods(mods_vf,mods_nonvar_vf):
    mod_ct_vf=helper_functions.normalize_counter(Counter(dict(mods_vf.value_counts())))
    plt.figure('discrepant peptide lengths')
    chist_vf=pd.DataFrame.from_dict(dict(mod_ct_vf.most_common(20)),orient='index')
    figname='nonvariant_unexpected_mods.png'
    title='All unexpected modifications found in non-variant peptides'
    mod_ct_non_var_vf=helper_functions.normalize_counter(Counter(dict(mods_nonvar_vf.value_counts())))
    chist_non_var_vf=pd.DataFrame.from_dict(dict(mod_ct_non_var_vf.most_common(20)),orient='index')
    combi=pd.concat([chist_non_var_vf,chist_vf],axis=1,sort=False).fillna(0)
    combi.columns=['All nonvariant VF', 'Mislabeled as nonvariant VF']
    combi=combi[~((combi['All nonvariant VF']<2)&(combi['Mislabeled as nonvariant VF']<2))]
    combi.plot(kind='bar',title=title)
    plt.ylabel("% Peptides")
    plt.xlabel("PTMs")
    plt.tight_layout()
    plt.savefig(figname)
    plt.clf()
    plt.close()

def plot_peplengths(lenct_vc,lenct_var,lenct_nonvar_vc,variant=False):
    labels=['Variant VC','All variant peptides','Non-variant peptides']
    figname='variant_peptide_length.png'
    plt.figure('discrepant peptide lengths')
    sns.set_context('paper')
    chist_vc=pd.DataFrame.from_dict(lenct_vc,orient='index').sort_index()
    chist_var=pd.DataFrame.from_dict(lenct_var,orient='index').sort_index()
    chist_nonvar_vc=pd.DataFrame.from_dict(lenct_nonvar_vc,orient='index').sort_index()
    combi=pd.concat([chist_vc,chist_var,chist_nonvar_vc],axis=1,sort=True).fillna(0)
    combi.columns=labels
    combi.head(30).plot(kind="bar",title="Length distributions of variant and non-variant peptides")
    #stats to look at the difference between the 2 columns
    t2, p2 = stats.ttest_ind(combi['Non-variant peptides'],combi['Variant VC']) #related or independent samples? (rel or ind?)
    print("peptide length statistics (independent t-test)")
    print("t = " + str(t2))
    print("p = " + str(p2))
    # combi=combi.reindex(new_index)
    plt.ylabel("Percentage")
    plt.xlabel("Length peptide")
    plt.legend(loc='upper right')
    plt.savefig(figname)
    plt.clf()
    plt.close()

def plot_final_venns(df_vc,df_vf):
    #get appropriate counters
    vc=Counter(dict(df_vc['variant_peptide'].value_counts()))
    vf=Counter(dict(df_vf['variant_peptide'].value_counts()))
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

