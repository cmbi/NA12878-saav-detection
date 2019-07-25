#!/usr/bin/env python3

'''this script will replace the protein field of the ionbot output with inferred proteins from percolator
will also concatenate all the output to a single file'''

import os, sys
import multiprocessing as mp
import pandas as pd


def main_replacement(percolatordir,csvname,ibcsvpath):
    # percolator_directory=os.fsencode(percolatordir)
    percolator_outfile=os.path.join(percolatordir,csvname)
    percolator_df=pd.read_csv(percolator_outfile,sep='\t')
    percolator_df.columns=["PSMId","proteinIds"]
    #take subset of the percolator file with that csv name in it
    temp=read_df_in_chunks(os.path.join(ibcsvpath,csvname),percolator_df,csvname, 1000)
    # temp["scan_id"]=temp["scan_id"].astype(str)+'_'+csvname #make scan ids unique again when concatenating all files
    print(csvname) #keep track of progress
    return(temp)

def read_df_in_chunks(directory,percolator_subset,csvname, chunksize):
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv(directory, chunksize=chunksize) # chunksize represents number of rows read per chunk
    chunk_list = []  # append each chunk df here 
    # Each chunk is in df format
    for chunk in df_chunk:  
        # perform data filtering/processsing
        chunk_filter = chunk_preprocessing(chunk,csvname)
        chunk_processed = replace_proteins(chunk_filter,percolator_subset)
        # Once the data filtering is done, append the chunk to list
        chunk_list.append(chunk_processed)
    # concat the list into dataframe 
    df_concat = pd.concat(chunk_list) # this is your final dataframe
    return(df_concat)

def chunk_preprocessing(df_chunk,csvname):
    new_chunk=df_chunk[(df_chunk['ri_126.1277']>0)]
    new_chunk=new_chunk[['scan_id','charge','precursor_mass','matched_peptide','modifications','ionbot_psm_score','DB','unexpected_modification','ms2pip_pearsonr','proteins','num_unique_pep_ids']]
    new_chunk["scan_id"]=new_chunk["scan_id"].astype(str)+'_'+csvname
    return(new_chunk)

def fetch_proteins(scanid,protinfdf):
    '''return the proteins as inferred by percolator algorithm "best-scoring peptide" instead of default ionbot no inference'''
    if len(protinfdf.loc[protinfdf['PSMId']==scanid])>0:
        sub_df=protinfdf.loc[protinfdf['PSMId']==scanid,'proteinIds']
        return(str(sub_df.iloc[0]))
    return('')

def replace_proteins(ibout_df, percolator_out_df):
    '''iterate through ib chunk df and replace the protein field with the appropriate protein field from the percolator output'''
    new_ib_df=pd.DataFrame(columns=['scan_id','charge','precursor_mass','matched_peptide','modifications','ionbot_psm_score','DB','unexpected_modification','ms2pip_pearsonr','proteins','num_unique_pep_ids'])
    for row in ibout_df.iterrows():
        scanid=row[1][0]
        prot_ids=fetch_proteins(scanid,percolator_out_df)
        if prot_ids!='':
            row[1][9]=prot_ids
        new_ib_df=new_ib_df.append(row[1])
    return(new_ib_df)

def file_list(path):
    directory= os.fsencode(path)
    list_files=[]
    for csvfile in os.listdir(directory):
        csvname=os.fsdecode(csvfile)
        if csvname.endswith('.csv'):
            list_files.append(csvname)
    return(list_files)

def main(ibcsvpath,percolatordir,outfile):
    '''do replacements per csv file to speed up the search, concatenates it all at the end'''
    pool=mp.Pool(processes=20)
    list_csvs=file_list(ibcsvpath)
    results=[pool.apply_async(main_replacement,args=(percolatordir,x,ibcsvpath)) for x in list_csvs]
    output = [p.get() for p in results]
    ionbotout=pd.concat(output)
    ionbotout.to_csv(outfile,index=False)
    return('finished')

#run twice, one for varfree and one for varcontaining
main(sys.argv[1],sys.argv[2],sys.argv[3])