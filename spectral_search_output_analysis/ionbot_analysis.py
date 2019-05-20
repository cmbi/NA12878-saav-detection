#!/usr/bin/env python3

import matplotlib, re, os
from matplotlib_venn import venn3, venn3_unweighted
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter

# Set the visualization settings
matplotlib.rcParams['axes.titlesize'] = 'xx-large'
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
from sets import Set

'''
this script will concatenate all the separate spectral match files generated from ionbot and process them to answer research questions

questions:
- horizontal/vertical coverage of the whole proteome?
- extent of multiple mapping
-
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

def import_cpdt(cpdt):
    cpdt_pep={}
    full_seqs={}
    with open(cpdt) as c:
        add=False
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                if key in set(proteins_covered):
                    add=True
                    cpdt_pep[key]={}
                    full_seqs[key]
                else:
                    add=False
            elif add and 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(':')[0]
                cpdt_pep[key][lp]=0
            elif add and 'PEPTIDE' not in line:
                full_seqs[key]=line.strip()
    return(cpdt_pep, full_seqs)



    
    
    
    
    
    
    
    
    
    
    