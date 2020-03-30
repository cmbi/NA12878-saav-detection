#!/usr/bin/env python3
import os, re
import pandas as pd
import helper_functions
import dask.dataframe as dd
# import glob


def concatenate_csvs(csvpath,vf=False):
    # directory= os.fsencode(csvpath)
    values={'unexpected_modification':'none'}
    ionbotout=dd.read_csv(f"{csvpath.strip('/')}/*.mgf.ionbot.csv")
    ionbotout=ionbotout[ionbotout['ri_126.1277']>0]
    # ionbotout['scan_id']=ionbotout['title'].str.split(' ').apply(lambda x: x[0]) #get unique identifier in the title column
    ionbotout=ionbotout.drop(['ri_126.1277','ri_127.1311','ri_128.1344','ri_129.1378','ri_130.1411','ri_131.1382'],axis=1)
    ionbotout["DB"]=ionbotout["DB"].map({'D':True,'T':False})
    ionbotout['peptide']=ionbotout['matched_peptide'].str.replace('I|L','x',regex=True)
    ionbotout['source_dict']=ionbotout['proteins'].apply(helper_functions.bin_hits_by_source,meta=str)
    if vf:
        ionbotout['pred_aa_sub']=ionbotout['modifications'].apply(lambda x: re.findall('[A-Z]{1}[a-z]{2}->[A-Z]{1}[a-z]{2}\[[A-Z]{1}\]',x),meta=str).apply(lambda y: re.findall('[A-Z]{1}[a-z]{2}',y[0]) if len(y)>0 else '',meta=str).apply(helper_functions.sub_conversion,meta=str)
    return(ionbotout.compute(scheduler='processes',num_workers=30))

def il_sensitive_read_csv(csvpath,to_replace=['peptide','ref_counterpart']):
    '''read in the variant peptides and counterparts, replace 
    '''
    #names may change if add substitution/variant status: ['protein','variant','counterpart','sub','start','chr','genomic_pos','is_het']
    df=pd.read_csv(csvpath)
    for col in to_replace:
        df[col]=df[col].replace(to_replace='I|L',value='x',regex=True)
    df['substitution']=df['substitution'].str.replace('I','X').str.replace('L','X')
    return(df.rename({'peptide':'variant_peptide'},axis=1))

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

# def import_cpdt(cpdt,fullSeq):
#     ''' read the cpdt files into a data structure
#     this function can also handle the cpdt files generated with interesting_peptide_finder (only peptides with SNVs)
#     {protein_ID:{pep1:0, pep2:0, pep3:0}} #counts
#     {pep1:prob, pep2:prob, pep3:prob} #probabilities
#     '''
#     cpdt_pep={}
#     cpdt_probs={}
#     full_seqs={}
#     with open(cpdt) as c:
#         for line in c:
#             if line.startswith('>'):
#                 key=line.strip()[1:]
#                 key=helper_functions.get_id(key)
#                 cpdt_pep[key]={}
#                 full_seqs[key]=''
#             elif 'PEPTIDE' in line:
#                 lp=line.split('PEPTIDE ')[1]
#                 lp=lp.split(':')
#                 cpdt_pep[key][lp[0]]=0
#                 cpdt_probs[lp[0]]=float(lp[1].strip())
#             elif 'PEPTIDE' not in line:
#                 full_seqs[key]=line.strip()
#     if fullSeq:
#         return(cpdt_pep, full_seqs)
#     return(cpdt_pep,cpdt_probs)

def import_gff(gfffile,isBed):
    '''use the gfffile to associate what proteins belong to which chromosome, in order to show the chromosomal distribution
    '''
    origininfo={'transcript_id':[],'chromosome':[],'strand':[]}
    with open(gfffile) as handle:
        for line in handle:
            info=line.split('\t')
            if len(info)>5:
                if not isBed: # if true gff
                    if 'transcript_type=protein_coding' in line:
                        tid=line.split('transcript_id=')[1]
                        tid=tid.split(';')[0]
                        origininfo['transcript_id'].append(tid)
                        origininfo['chromosome'].append(info[0])
                        origininfo['strand'].append(info[6])
                else:
                    origininfo['transcript_id'].append(info[3])
                    origininfo['chromosome'].append(info[0])
                    if info[5]!='+' and info[5]!='-':
                        origininfo['strand'].append('unknown')
                    else:
                        origininfo['strand'].append(info[5])
    return(pd.DataFrame(data=origininfo))


def create_chromosome_reference(gfffile,bedfile):
    info_ref=import_gff(gfffile,False) #import gff3 file annotations from gencode
    info_ont=import_gff(bedfile,True) #import bed file annotations from ont (converted from psl)
    return(pd.concat([info_ref,info_ont]))

# def import_cpdt_simple(cpdt):
#     cpdt_pep={}
#     with open(cpdt) as c:
#         for line in c:
#             if line.startswith('>'):
#                 key=line.strip()[1:]
#                 key=helper_functions.get_id(key)
#                 cpdt_pep[key]=[]
#             elif 'PEPTIDE' in line:
#                 lp=line.split('PEPTIDE ')[1]
#                 lp=lp.split(':')[0]
#                 cpdt_pep[key].append(lp)
#     return(cpdt_pep)

