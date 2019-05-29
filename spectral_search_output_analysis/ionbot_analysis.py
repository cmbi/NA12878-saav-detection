#!/usr/bin/env python3

import matplotlib, re, os
from matplotlib_venn import venn3, venn3_unweighted
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sets import Set
from collections import Counter

# Set the visualization settings
matplotlib.rcParams['axes.titlesize'] = 'xx-large'
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

'''
this script will concatenate all the separate spectral match files generated from ionbot and process them to answer research questions

questions:
- horizontal/vertical coverage of the whole proteome?
- extent of multiple mapping, proportion of mapping that is exclusive to the long-read transcriptome
- are my peptides of interest present? How often does ionbot correctly predict an SNV?
'''
def concatenate_csvs(folder_path):
    directory= os.fsencode(csvpath)
    ionbotout=pd.DataFrame()
    for csvfile in os.listdir(directory):
        csvname=os.fsdecode(csvfile)
        if csvname.endswith('.csv'):
            temp=read_df_in_chunks(csvname, 1000)
            ionbotout=pd.concat([ionbotout,temp])
    return(ionbotout)

def read_df_in_chunks(filename, chunksize):
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv(filename, chunksize=chunksize) # chunksize represents number of rows read per chunk

    chunk_list = []  # append each chunk df here 

    # Each chunk is in df format
    for chunk in df_chunk:  
        # perform data filtering 
        chunk_filter = chunk_preprocessing(chunk)

        # Once the data filtering is done, append the chunk to list
        chunk_list.append(chunk_filter)

    # concat the list into dataframe 
    df_concat = pd.concat(chunk_list) # this is your final dataframe
    
    return(df_concat)

def chunk_preprocessing(df_chunk):
    new_chunk=df_chunk[(df_chunk['ri_126.1277']>0) & (df_chunk['q_value']<=0.01) & (df_chunk['DB']=='T')]
    new_chunk=new_chunk[['scan_id','charge','precursor_mass','matched_peptide','modifications','ionbot_psm_score','DB','unexpected_modification','ms2pip_pearsonr','proteins','num_unique_pep_ids']]
    return(new_chunk)

def import_coding_transcriptids(sources):
    '''
    input: paths of protein sequences from reference and cell-specific transcriptome translation
    output: ids of the combination
    
    '''
    transcript_ids=[]
    sources=['/data/data/genome_grch38/gencode.v29.pc_translations.fa','/data/data/spectra/dictionary/dumb_orfprediction_setA/flair.setA.final.pep']
    for f in sources:
        with open(f) as handle:
            for line in handle:
                if line.startswith('>'):
                    if ' ' in line.strip():
                        tid=line.split(' ')[0]
                        transcript_ids.append(tid[1:])
                    else:
                        transcript_ids.append(line.strip()[1:])
    return(transcript_ids)

def detected_proteins(row,pco):
    proteins_covered=pco
    if '||' in row[1][9]:
        ids=row[1][9].split('||')
        for idu in ids:
            proteins_covered[idu]+=1
    else:
        proteins_covered[row[1][9]]+=1
    return(proteins_covered)

def import_cpdt(cpdt,wantFull):
    ''' read the cpdt files into a data structure
    this function can also handle the cpdt files generated with interesting_peptide_finder (only peptides with SNVs)
    {protein_ID:{pep1:0, pep2:0, pep3:0}}
    '''
    cpdt_pep={}
    full_seqs={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                if '|m.' in key:
                    key=key.split('|m.')[0]
                cpdt_pep[key]={}
                full_seqs[key]=''
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                if ':' in lp:
                    lp=lp.split(':')[0]
                cpdt_pep[key][lp]=0
            elif 'PEPTIDE' not in line:
                full_seqs[key]=line.strip()
    if wantFull:
        return(cpdt_pep, full_seqs)
    return(cpdt_pep)

def import_gff(gfffile,isBed):
    '''use the gfffile to associate what proteins belong to which chromosome, in order to show the chromosomal distribution
    '''
    chromdict={}
    with open(gfffile) as handle:
        for line in handle:
            info=line.split('\t')
            if len(info)>5:
                if not isBed:
                    if 'transcript_type=protein_coding' in line:
                        tid=line.split('transcript_id=')[1]
                        tid=tid.split(';')[0]
                        chromdict[tid]=info[0]
                else:
                    chromdict[info[3]]=info[0]
    return(chromdict)

def match_peps(row,cpdt_pep,olddetected):
    '''match predicted mutated peptides to observed
    how many unique variant peptides are detected, from how many unique proteins?
    how many total instances of correct/incorrect variant peptides are detected?
    build up the dictionary with every iteration of the 
    '''
    prot=row[1][9]
    pep=row[1][3]
    mod=str(row[1][7])
    detected=olddetected
    if '|m.' in prot:
        prot=prot.split('|m.')[0]
    if prot in cpdt_pep:
        for p,ct in cpdt_pep.items():
            if p==pep:
                if prot not in detected:
                    detected[prot]={}
                detected[prot][pep]=ct
    return(detected)

def coverage_measure(cpdt_pep):
    high_cov_vert={}
    high_cov_hor={}
    perc_cov_dist=[]
    vert_cov=[]
    for p,peps in cpdt_pep.items():
        seq=full_seqs[p]
        remains=seq
        count_pep=0
        for s,c in peps.items():
            if c>5: #what constitutes a "true" hit
                count_pep+=c
                if s in remains:
                    remains=remains.replace(s,'')
                else:
                    prefix=re.split('R|K',s)
                    for p in prefix:
                        if len(p)>3 and p in remains:
                            remains=remains.replace(p,'')
        perc_cov=float((len(seq)-len(remains))/len(seq))*100
        perc_cov_dist.append(perc_cov)
        vert_cov.append(count_pep)
        if perc_cov>50:
            high_cov_hor[p]=peps
        if count_pep>100:
            high_cov_vert[p]=peps
    return(high_cov_hor,high_cov_vert,perc_cov_dist)

def bin_hits_by_source(row,oro,ono,ob):
    '''sort peptide hits by their source dictionary'''
    ref_only=oro
    ont_only=ono
    both=ob
    if '||' in row[1][9]:
        ids=row[1][9].split('||')
        ont=False
        ref=False
        for i in ids:
            if '|m.' in i: #assumes that the proteins from the reference were not generated with ANGEL
                ont=True
            else:
                ref=True
        if ont and ref:
            both.add(row[1][0])
        elif ont:
            ont_only.add(row[1][0])
        elif ref:
            ref_only.add(row[1][0])
        else:
            raise Exception('Unexpected protein found')
    else:
        if 'HUMAN' in (row[1][9]):
            ref_only.add(row[1][0])
        else:
            ont_only.add(row[1][0])
    return(ref_only,ont_only,both)

def plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_combi):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    plt.figure("Pearson R distribution")
    matplotlib.rcParams['axes.titlesize'] = 'xx-large'
    matplotlib.rcParams['axes.labelsize'] = 'x-large'
    matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
    sns.distplot(ibdf_refonly['ms2pip_pearsonr'], hist=False, label='Human reference only',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_ontonly['ms2pip_pearsonr'], hist=False, label='Transcriptome translation only',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_combi['ms2pip_pearsonr'], hist=False, label='Reference + transcriptome translation',axlabel='Pearson R Correlation')
    plt.legend(prop={'size':30})
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    return("Scores plot made")

def plot_source_piechart():
    '''this function will plot the source piechart and save it to a pdf'''
    return 0

def plot_coverage_plots(cpdt_pep):
    '''this function will plot the graphs that correspond to the coverage of the proteome
    - vertical coverage
    - horizontal coverage
    - chromosome distribution
    '''
    return 0

def make_report(hits_df):
    '''this function will output a report with general information collected about the hits
    which includes:
    - how many unique peptides are found
    - how many unique proteins are found (and how many are in transcriptome v reference)
    - inventory PTMs (including how many SNVs there are)
    '''
    uniquepeptides=hits_df['peptide'].nunique()
    return 0

    
def main(directory_ontonly, directory_refonly, directory_combination):
    '''this is the main function that will iterate over the giant pandas df and perform all analyses and make all figures
    this is written in a way that iteration should only be done once.
    '''
    #imports
    ibdf_ontonly=concatenate_csvs(directory_ontonly)
    ibdf_refonly=concatenate_csvs(directory_refonly)
    ibdf_combi=concatenate_csvs(directory_combination)
    cpdt_pep,full_seqs=import_cpdt(cpdtfile,True) #import cpdt will all peptides (cat gencode and flair beforehand)
    mut_cpdt_pep=import_cpdt(cpdtfile_mut,False) #import the cpdt file with all mutant peptides
    chromdict_ref=import_gff(gfffile,False) #import gff3 file annotations from gencode
    chromdict_ont=import_gff(bedfile,True) #import bed file annotations from ont (converted from psl)
    chromdict={**chromdict_ont,**chromdict_ref} #combine the 2 dictionaries

    
    #initialize data structures to collect information
    proteins_covered=Counter() #proteins detected
    detected={} #variant peptides detected
    ref_only=set() #scan ids in the reference set
    ont_only=set() #scan ids in the ont set
    both=set() #scan ids that matched to both ref and ont proteins

    #qc function


    #iterate to fill the data structures
    for row in ibdf.iterrows():
        proteins_covered=detected_proteins(row,proteins_covered) #what proteins from the proteome are covered and in what amounts
        ref_only,ont_only,both=bin_hits_by_source(row,ref_only,ont_only,both) #what dictionaries do the hits come from
        detected=match_peps(row,mut_cpdt_pep,detected) #what mutant peptides are detected, how many, and what proteins they come from

    #create
    return 0
    
    
    
    
    
    
    
    