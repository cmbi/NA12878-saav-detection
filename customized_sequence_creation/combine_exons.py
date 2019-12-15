#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import numpy as np
#from enum import Enum
from tqdm import tqdm
import gffpandas.gffpandas as gffpd
import multiprocessing as mp
import argparse
import traceback

'''
combine exons: given "phased" exon sequences from phase_exons, assemble full length transcripts (per haplotype)
overview functions:
    -read_gff3_into_frame: this will take either a gff3 or bed file and put it into a pandas dataframe with the relevant information for this analysis
    -child_initialize: globalizes and initializes the exon dataframe for processing
    -worker_process: the main function to assemble a transcript, run once per transcript
    -variant_renumber: renumber variants based on strandedness and length of existing sequence
    -reverse_complement: in case minus strand, make reverse complement
        -make_complement: does the reverse-complimenting
    -filter_nones: removes the "None"s in a list of variants

in:
    -number of cpu to use
    -gff3 or bed file of transcript junctions
    -output file of "phase_exons" - tab-seperated file of exon information and sequences
    -optional- cds file in gff3 format, in the case when cds annotated with genomic coordinates
    -names of output files

returns: 
    -a fasta file of assembled transcripts that include variants and are phased (variant-free transcripts not written)
    -a tab seperated report of what variants are in each transcript
        columns: transcript, haplotype, het variant positions, het variant origins, hom variant positions, hom variant origins, start codon positions

'''
                    
def read_gff3_into_frame(gff):
    '''helper function to combine_exons
    input: junction file (gff), and T/F depending on whether you are actually putting in a 0-based (bed) file. True for bed, False for gff.
    this function assumes that the bed file chr is in the foormat "chr2" while gff is in format "2"
    '''
    if gff[-3:]=='bed':
        df=pd.read_csv(gff,sep='\t',header=None)
        cols_to_drop=[3,4,6,7,8]
        transcript_seperator='transcript:'
        exon_position_seperator='rank='
        shift=0
    elif gff[-4:]=='gff3':
        df=gffpd.read_gff3(gff).df
        cols_to_drop=[1,2,5,7]
        transcript_seperator='transcript_id='
        exon_position_seperator='exon_number='
        shift=-1
    else:
        raise Exception("Wrong junction file type input - must be gff3 or bed file")
    df.drop(df.columns[cols_to_drop],axis=1,inplace=True) #maybe these shouldn't be hardcoded...
    df.columns=['chromosome','start','end','sense','transcript']
    df['start']=df['start'].astype(int)+shift
    df['transcript']=df['transcript'].str.split(transcript_seperator).apply(lambda x: x[1])
    df['transcript'],df['rank']=df['transcript'].str.split(exon_position_seperator,1).str
    df['transcript']=df['transcript'].str.split(';',1).apply(lambda x: x[0])
    df['rank']=df['rank'].str.split(';',1).apply(lambda x: x[0])
    if 'start_codon' in gff: #FOR START CODONS 
        df.drop(columns=['sense','chromosome'],inplace=True) #drop redundant
        df=df.groupby(['transcript','rank']).agg(lambda x:list(x))#.agg(dict(start=lambda x: ','.join(str(i) for i in x))) #get list of start codons in case more than one
        df.rename({'start':'start_codon_start','end':'start_codon_end'},axis=1,inplace=True) #rename start codon
    return(df)

def reverse_complement(df):
    '''read in exons and process the minus strand sequences
    chromosome      start   end     haplotype       pos_het het_org pos_hom hom_org sequence
    '''
    reverse_strand=df[df['sense']=='-']
    forward_strand=df[df['sense']=='+']
    reverse_strand['sequence']=reverse_strand['sequence'].apply(lambda x: make_complement(str(x)))
    # reverse_strand['start_codon_shift']=reverse_strand.apply(lambda x: start_codon_coord(x['start_codon_end'],x['end'],reverse=True),axis=1)
    # forward_strand['start_codon_shift']=forward_strand.apply(lambda x: start_codon_coord(x['start_codon_start'],x['start']),axis=1)
    return(pd.concat([forward_strand,reverse_strand],axis=0, ignore_index=True))

def make_complement(sequence):
    '''helper function to combine_exons
    makes complimenary mRNA in case of antisense strand'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
    return(reverse_complement)

def filter_nones(the_list):
    '''Remove "None"s within the list of variants, which sometimes happens if exons with no variants are added between exons that do have variants'''
    new_list=[value for value in the_list if value != "None"]
    if not new_list:
       return(str(None))
    return(new_list)

def start_codon_coord(start_codon_pos_list,exon_pos,reverse=False):
    '''renumber the start codon based on the exon start
    note: for reverse strand, you should input the "end" instead of "start"'''
    if reverse:
        return([int(exon_pos) - int(x) for x in start_codon_pos_list])
    return([int(x) - int(exon_pos) for x in start_codon_pos_list])


def variant_renumber(com_seperated_variantlist,len_whole_seq,strand,len_exon):
    '''takes a comma seperated string of variants and the length of the whole sequence
    renumbers the variants to adjust for the minus strand and adds the length of the existing sequence
    input: comma-seperated list of variants e.g. "23,634,2335", length of existing seq, strandedness and length of exon
    output: new list of variants in the same format, but renumbered
    '''
    if com_seperated_variantlist!="None":
        assert type(com_seperated_variantlist)==str, com_seperated_variantlist
        if ',' in com_seperated_variantlist:
            com_seperated_variantlist=com_seperated_variantlist.split(',')
            new_var_pos=[]
            for variant in com_seperated_variantlist:
                if strand=='-':
                    variant=str(len_exon-int(variant))
                new_var_pos.append(str(len_whole_seq+int(variant)))
            return(','.join(new_var_pos))
        else:
            if strand=='-':
                com_seperated_variantlist=str(len_exon-int(com_seperated_variantlist))
            return(str(len_whole_seq+int(com_seperated_variantlist)))
    return(com_seperated_variantlist)

def worker_process(transcript):
    ''' The main task for assembling transcripts
    Assembles transcripts using exon order according to the junction file
    Renumbers variant positions according to position in the total transcript sequence

    Input: transcript name
    output: lines to print to the output fasta/report
    '''
    #store the sequence as it's being built
    sequence_zero,sequence_one='',''
    #booleans to check for variants
    hasHeterozygous,hasHomozygous=False,False
    #homologous and heterozgyous variants, and their original position in the genome. Record 0 and 1 haplotype seperately in tuple
    var_hom,hom_org,var_het,het_org=([],[]),([],[]),([],[]),([],[])
    #keep track of original sequence length
    len_seq_org=0
    try:
        trans_info=exon_df.loc[exon_df["transcript"]==transcript]#get exons in transcript
        #check for consecutive
        int_ranks=sorted([int(x) for x in list(trans_info['rank'].unique())])
        if int_ranks != list(range(min(int_ranks), max(int_ranks)+1)): 
            print(int_ranks)
            raise Exception(f"ranks for transcript {transcript} are inconsistent; exon may be missing")
        #get the start pos
        # fiveprime=trans_info.loc[trans_info["rank"]=='1','five_prime_shift'].iloc[0]
        for rank in int_ranks: #iterate through each exon by its rank, by ascending
            rank=str(rank)
            seq_rows=trans_info[trans_info["rank"]==rank].reset_index(drop=True)
            if len(seq_rows.index)==1: #if no heterozygous variants in exon
                assert (seq_rows.iloc[0]['pos_het']=="None"), f"there are heterozygous positions at {seq_rows.iloc[0]['pos_het']}"
                if seq_rows.iloc[0]["pos_hom"]!="None":
                    hasHomozygous=True
                var_hom[0].append(variant_renumber(seq_rows.iloc[0]["pos_hom"],len(sequence_zero),seq_rows.iloc[0]['sense'],len(seq_rows.iloc[0]['sequence'])))
                var_hom[1].append(variant_renumber(seq_rows.iloc[0]["pos_hom"],len(sequence_one),seq_rows.iloc[0]['sense'],len(seq_rows.iloc[0]['sequence'])))
                hom_org[0].append(seq_rows.iloc[0]["hom_org"]) #transfer over the origins
                hom_org[1].append(seq_rows.iloc[0]["hom_org"])
                #finally, add the exon sequence to the transcript
                sequence_zero+=seq_rows.iloc[0]["sequence"]
                sequence_one+=seq_rows.iloc[0]["sequence"]
            elif len(seq_rows.index)==2: #if heterozygous variants in exon
                hasHeterozygous=True
                #separate the 2 exon haplotypes
                hap_zero_exon_row=seq_rows[seq_rows["haplotype"]=="0"]
                hap_one_exon_row=seq_rows[seq_rows["haplotype"]=="1"]
                #collect the corresponding variants
                if hap_zero_exon_row.shape[0]!=1 or hap_one_exon_row.shape[0]!=1:
                    return('',[]) #this happened in 2 cases
                if hap_zero_exon_row.iloc[0]["pos_het"]=="None" or hap_one_exon_row.iloc[0]["pos_het"]=="None":
                    print(seq_rows)
                    raise Exception("Missing heterozygous variant positions in a heterzygous exon")
                #collect all information
                var_het[0].append(variant_renumber(hap_zero_exon_row.iloc[0]["pos_het"],len(sequence_zero),hap_zero_exon_row.iloc[0]['sense'], \
                len(hap_zero_exon_row.iloc[0]['sequence'])))
                var_het[1].append(variant_renumber(hap_one_exon_row.iloc[0]["pos_het"],len(sequence_one),hap_one_exon_row.iloc[0]['sense'], \
                len(hap_one_exon_row.iloc[0]['sequence'])))
                het_org[0].append(hap_zero_exon_row.iloc[0]["het_org"])
                het_org[1].append(hap_one_exon_row.iloc[0]["het_org"])
                var_hom[0].append(variant_renumber(hap_zero_exon_row.iloc[0]["pos_hom"],len(sequence_zero),hap_zero_exon_row.iloc[0]['sense'], \
                len(hap_zero_exon_row.iloc[0]['sequence'])))
                var_hom[1].append(variant_renumber(hap_one_exon_row.iloc[0]["pos_hom"],len(sequence_one),hap_one_exon_row.iloc[0]['sense'], \
                len(hap_one_exon_row.iloc[0]['sequence'])))
                hom_org[0].append(hap_zero_exon_row.iloc[0]["hom_org"])
                hom_org[1].append(hap_one_exon_row.iloc[0]["hom_org"])
                #finally, add the exon sequences to the transcript
                sequence_zero+=hap_zero_exon_row.iloc[0]["sequence"]
                sequence_one+=hap_one_exon_row.iloc[0]["sequence"]
            else:
                print(trans_info)
                raise Exception('exon not found')
        #process and return full sequences
        if hasHeterozygous: #if there was a heterozygous variant detected somewhere in the transcript
            return(f">{transcript}_h0\n{sequence_zero}\n>{transcript}_h1\n{sequence_one}\n", \
            [f"{transcript}\t0\t{','.join(var_het[0])}\t{','.join(het_org[0])}\t \
            {','.join(filter_nones(var_hom[0])) if hasHomozygous else str(None)}\t \
            {','.join(filter_nones(hom_org[0])) if hasHomozygous else str(None)}\n",
            f"{transcript}\t1\t{','.join(var_het[1])}\t{','.join(het_org[1])}\t \
            {','.join(filter_nones(var_hom[1])) if hasHomozygous else str(None)}\t \
            {','.join(filter_nones(hom_org[1])) if hasHomozygous else str(None)}\n"])
        elif hasHomozygous:
            return(f">{transcript}\n{sequence_zero}\n", \
            [f"{transcript}\tNone\t{str(None)}\t{str(None)}\t{','.join(filter_nones(var_hom[0]))}\t \
            {','.join(filter_nones(hom_org[0]))}\n"])
        # else: #no variants= don't return! we only want the sequences that have variants, otherwise creating duplicates
        #     return()
    except Exception as e:
        raise e

def print_results(results,outfasta,report):
    '''
    takes the exons created from the phase_exons function and stitches them together according to the gff fasta to make full transcripts
    
    ''' 
    f=open(outfasta,'w')
    r=open(report,'w')
    output = [p.get() for p in results] #list of tuplies
    #write transcript to fasta
    r.writelines('\t'.join(['transcript','haplotype','pos_het','het_origin','pos_hom','hom_origin'])+'\n')
    for o in output:
        if o:
            f.write(o[0])
            for line in o[1]:
                r.write(line)
    r.close()
    f.close()
    return('done')

def child_initialize(_exon_df):
     global exon_df
     exon_df = _exon_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--cpu',help='CPU number',required=True)
    parser.add_argument('--exons',help='Phased and variant replaced exons, tab seperated from phase_exons.py',required=True)
    parser.add_argument('--jun', help='file transcript junctions (GFF or BED accepted)', required=True)
    parser.add_argument('--out', help='Output fasta', required=True)
    parser.add_argument('--report', help='Output file report', required=True) #in favor of reporting the variants in the header of the fasta file
    parser.add_argument('--debug', help='debug option',action='store_true', required=False) #in favor of reporting the variants in the header of the fasta file
    parser.add_argument('--sc',help='start codons (if only available in the genomic coordinates)',required=False)
    args=vars(parser.parse_args())
    print('Reading in haplotype-specific exons...')
    exon_frame=pd.read_csv(args['exons'],sep='\t')
    print('Reading in transcript annotations...')
    junction_frame=read_gff3_into_frame(args['jun']) #can put bed file 
    print("Preparing for analysis...")
    exon_df=reverse_complement(pd.merge(junction_frame, exon_frame, on=['chromosome','start','end']))
    # exon_df.replace({'None':np.nan}, inplace=True)
    transcripts=exon_df.transcript.unique()
    print(f'building {len(transcripts)} transcripts')
    if args['debug']:
        transcripts=['ENST00000337907.7']
        pool=mp.Pool(processes=1, initializer=child_initialize, initargs=(exon_df,))
    else:
        pool=mp.Pool(processes=int(args['cpu']), initializer=child_initialize, initargs=(exon_df,))
    results=[pool.apply_async(worker_process,args=(transcript,)) for transcript in transcripts] #each transcript processed with seperate cpu
    print_results(results,args['out'],args['report'])
