#!/usr/bin/env python3

#########
#idea: move data into SQL database and move the quality control plots above to another script to cut down on processing time, since the quality control steps are the only ones that need all, unfiltered data
#insert here an SQL query here to import data that meets preprocessing cutoffs
#would need to concatenate all csv files, change the scan ids to include the file names, and upload them onto an SQL server

def concatenate_csvs(csvpath):
    directory= os.fsencode(csvpath)
    ionbotout=pd.DataFrame()
    for csvfile in os.listdir(directory):
        csvname=os.fsdecode(csvfile)
        if csvname.endswith('.csv'):
            temp=read_df_in_chunks(os.path.join(csvpath,csvname), 1000)
            temp["scan_id"]=temp["scan_id"].astype(str)+'_'+csvname #make scan ids unique again when concatenating all files
            temp["title"]=temp["title"].astype(str).split('""')[1]
            ionbotout=pd.concat([ionbotout,temp])
    return(ionbotout)

def read_df_in_chunks(directory, chunksize):
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv(directory, chunksize=chunksize) # chunksize represents number of rows read per chunk
    chunk_list = []  # append each chunk df here 
    # Each chunk is in df format
    for chunk in df_chunk:  
        # perform data filtering 
        #chunk_filter = chunk_preprocessing(chunk)
        # Once the data filtering is done, append the chunk to list
        chunk_list.append(chunk)
    # concat the list into dataframe 
    df_concat = pd.concat(chunk_list) # this is your final dataframe
    return(df_concat)

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

    # Score distribution plot
    df['DB']=df['DB']=='D'
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
    plt.xlim(0.00001, 1)
    plt.legend()
    plt.title('q-value comparison search dictionaries')
    plt.tight_layout()
    plt.savefig('qval_comparison.png')


def plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_combi):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution")
    sns.distplot(ibdf_refonly['percolator_psm_score'], hist=False, label='Human reference only',axlabel='Percolator score')
    sns.distplot(ibdf_ontonly['percolator_psm_score'], hist=False, label='Transcriptome translation only',axlabel='Percolator score')
    sns.distplot(ibdf_combi['percolator_psm_score'], hist=False, label='Reference + transcriptome translation',axlabel='Percolator score')
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
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='T']['percolator_psm_score'], hist=False, label='Reference + transcriptome translation (no variants)',axlabel='Percolator score')
    sns.distplot(ibdf_combi_pg[ibdf_combi_pg['DB']=='T']['percolator_psm_score'], hist=False, label='Reference + transcriptome translation (with variants)',axlabel='Percolator score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr_pgvsopenmut.png")
    plt.close()
    return("Scores plot made")

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
    plots.plot_target_decoy(ibdf_combi.dropna(),"qc_pearsonr_decoy_varfree.png", plot_title="Search result variant-free")
    plots.plot_target_decoy(ibdf_combi_pg.dropna(),"qc_pearsonr_decoy_varcont.png", plot_title="Search result variant-containing")
    plots.plot_qvalues_comparison({'ONT only':ibdf_ontonly,'Ref only':ibdf_refonly,'Combi variant-containing':ibdf_combi_pg,'Combi variant-free':ibdf_combi},fdr_levels=[0.01])
