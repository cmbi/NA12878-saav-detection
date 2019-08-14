#!/usr/bin/env python3
import os
import pandas as pd
import helper_functions


def concatenate_csvs(csvpath):
    directory= os.fsencode(csvpath)
    ionbotout=pd.DataFrame()
    for csvfile in os.listdir(directory):
        csvname=os.fsdecode(csvfile)
        if csvname.endswith('.csv'):
            temp=read_df_in_chunks(os.path.join(csvpath,csvname), 1000)
            temp["scan_id"]=temp["scan_id"].astype(str)+'_'+csvname.split('.')[0] #make scan ids unique again when concatenating all files
            temp["title"]=csvname.split('.')[0]
            temp["DB"]=temp["DB"].map({'D':True,'T':False})
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

def chunk_preprocessing(df_chunk):
    new_chunk=df_chunk[(df_chunk['ri_126.1277']>0) & (df_chunk['q_value']<=0.01) & (df_chunk['DB']==False)]
    new_chunk=new_chunk[['scan_id','charge','precursor_mass','matched_peptide','modifications','percolator_psm_score','DB','unexpected_modification','ms2pip_pearsonr','proteins','num_unique_pep_ids']]
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

def import_cpdt(cpdt,fullSeq):
    ''' read the cpdt files into a data structure
    this function can also handle the cpdt files generated with interesting_peptide_finder (only peptides with SNVs)
    {protein_ID:{pep1:0, pep2:0, pep3:0}} #counts
    {pep1:prob, pep2:prob, pep3:prob} #probabilities
    '''
    cpdt_pep={}
    cpdt_probs={}
    full_seqs={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=helper_functions.get_id(key)
                cpdt_pep[key]={}
                full_seqs[key]=''
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(':')
                cpdt_pep[key][lp[0]]=0
                cpdt_probs[lp[0]]=float(lp[1].strip())
            elif 'PEPTIDE' not in line:
                full_seqs[key]=line.strip()
    if fullSeq:
        return(cpdt_pep, full_seqs)
    return(cpdt_pep,cpdt_probs)

def import_gff(gfffile,isBed):
    '''use the gfffile to associate what proteins belong to which chromosome, in order to show the chromosomal distribution
    '''
    chromdict={}
    stranddict={}
    with open(gfffile) as handle:
        for line in handle:
            info=line.split('\t')
            if len(info)>5:
                if not isBed:
                    if 'transcript_type=protein_coding' in line:
                        tid=line.split('transcript_id=')[1]
                        tid=tid.split(';')[0]
                        chromdict[tid]=info[0]
                        stranddict[tid]=info[6]
                else:
                    chromdict[info[3]]=info[0]
                    if info[5]!='+' and info[5]!='-':
                        stranddict[info[3]]='unknown'
                    else:
                        stranddict[info[3]]=info[5]
    return(chromdict,stranddict)


def create_chromosome_reference(gfffile,bedfile):
    chromdict_ref,stranddict_ref=import_gff(gfffile,False) #import gff3 file annotations from gencode
    chromdict_ont,stranddict_ont=import_gff(bedfile,True) #import bed file annotations from ont (converted from psl)
    chromdict={**chromdict_ont,**chromdict_ref} #combine the 2 dictionaries
    stranddict={**stranddict_ont,**stranddict_ref} #combine the 2 dictionaries
    return(chromdict,stranddict)

def import_cpdt_simple(cpdt):
    cpdt_pep={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=helper_functions.get_id(key)
                cpdt_pep[key]=[]
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(':')[0]
                cpdt_pep[key].append(lp)
    return(cpdt_pep)

