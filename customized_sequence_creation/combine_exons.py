#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import numpy as np
#from enum import Enum
from Bio.Seq import Seq
from tqdm import tqdm
import multiprocessing as mp
import argparse

'''
combine exons: given "phased" exon sequences from phase_exons, put the exons back together (per haplotype) to create full length transcripts
    -read_exons_into_frame: parses the exon file and organizes them into a pandas dataframe
    -read_gff_into_frame: parses the junction files and organizes them into a pandas dataframe
    -make_complement: in case minus strand, make reverse complement
    -process_variant: separates the exon coordinate from the genome coordinate

returns: a fasta file of completed transcripts that include variants and phased. 
the fasta header of the output includes a report of what heterozygous/homozygous
variants are included in the transcript, what the coordinate is in the transcript,
and what the genome coordinate is

>TRANSCRIPT_NAME_h0 pos_het:66|chr1|4385634,79|chr1|4385764 pos_hom:122|chr1|4387777
{sequence}
'''
                        
def print_results(results,outfile):
    '''
    takes the exons created from the phase_exons function and stitches them together according to the gff file to make full transcripts
    
    ''' 
    f=open(outfile,'w')
    output = [p.get() for p in results] #list of tuplies
    #write transcript to fasta
    for o in output:
        f.write(o)
    return('done')

def worker_process(transcript):
    sequence_zero=''
    sequence_one=''
    to_add=''
    isHetero=False
    var_hom,hom_org,var_het,het_org=([],[]),([],[]),([],[]),([],[])
    output_fasta=[]
    output_report=[]
    trans_info=junction_frame[junction_frame["transcript"]==transcript] #get junction information for the transcript
    trans_info=trans_info.sort_values(by="rank") #sort them in the correct order
    for tidx,trow in trans_info.iterrows(): #iterate through each junction and get the sequence
        seq_rows=exon_frame[(exon_frame["chromosome"]==trow["chromosome"]) & (exon_frame["end"]==trow["end"]) & (exon_frame["start"]==trow["start"])]
        if len(seq_rows.index)==1: #if no heterozygous variants in exon
            assert(seq_rows.iloc[0]['pos_het'] == None)
            var_hom[0].append(variant_renumber(seq_rows.iloc[0]["pos_hom"],len(sequence_zero),strand_correction=False))
            var_hom[1].append(variant_renumber(seq_rows.iloc[0]["pos_hom"],len(sequence_one),strand_correction=False))
            # if v:
            #     var_hom[0].extend([str(int(x[0])+len_existingz)+x[1] for x in v])
            #     var_hom[1].extend([str(int(x[0])+len_existingo)+x[1] for x in v])
            #finally, add the exon sequence to the transcript
            sequence_zero+=seq_rows.iloc[0]["sequence"]
            sequence_one+=seq_rows.iloc[0]["sequence"]
        elif len(seq_rows.index)==2: #if heterozygous variants in exon
            isHetero=True
            #separate the 2 exon haplotypes
            seq_row_hz=seq_rows[seq_rows["haplotype"]==0]
            seq_row_ho=seq_rows[seq_rows["haplotype"]==1]
            #collect the corresponding variants
            assert(seq_row_hz.iloc[0]["pos_het"]!=None)
            assert(seq_row_ho.iloc[0]["pos_het"]!=None)
            #fetch sequence, renumber if minus strand, add exon sequence to transcript
            var_het[0].append(variant_renumber(seq_row_hz.iloc[0]["pos_het"],len(sequence_zero),strand_correction=False))
            var_het[1].append(variant_renumber(seq_row_ho.iloc[0]["pos_het"],len(sequence_one),strand_correction=False))
            var_hom[0].append(variant_renumber(seq_row_hz.iloc[0]["pos_hom"],len(sequence_zero),strand_correction=False))
            var_hom[1].append(variant_renumber(seq_row_ho.iloc[0]["pos_hom"],len(sequence_one),strand_correction=False))
            #finally, add the exon sequences to the transcript
            sequence_zero+=seq_row_hz.iloc[0]["sequence"]
            sequence_one+=seq_row_ho.iloc[0]["sequence"]
            #update and store variant info
            # var_het[0].extend([str(int(x[0])+len_existingz)+x[1] for x in vhz])
            # var_het[1].extend([str(int(x[0])+len_existingo)+x[1] for x in vho])
            # if vhomz:
            #     var_hom[0].extend([str(int(x[0])+len_existingz)+x[1]  for x in vhomz])
            # if vhomo:
            #     var_hom[1].extend([str(int(x[0])+len_existingo)+x[1] for x in vhomo])
        else:# can have another elif for 0, and then an 
            print(trans_info)
            raise Exception('exon not found')
    if isHetero: #if there was a heterozygous variant detected somewhere in the transcript
        for_fasta=f">{transcript}_h0\n{sequence_zero}\n>{transcript}_h1\n{sequence_one}\n"
        for_report=f""
        header_zero='>'+transcript+'_h0 pos_het:'+','.join(var_het[0]) # look into f strings- more pythonic
        header_one='>'+transcript+'_h1 pos_het:'+','.join(var_het[1])
        if len(var_hom[0])>0:
            header_zero+=' pos_hom:'+','.join(var_hom[0])
        if len(var_hom[1])>0:
            header_one+=' pos_hom:'+','.join(var_hom[1])
        l1=header_zero+'\n'+sequence_zero+'\n'
        l2=header_one+'\n'+sequence_one+'\n'
        output_fasta=[l1,l2]
    else: #if no heterozygous variant detected, then both 0 and 1 haplotypes are the same
        header_zero='>'+transcript
        if len(var_hom[0])>0:
            header_zero+=' pos_hom:'+','.join(var_hom[0])
        l=header_zero+'\n'+sequence_zero+'\n'
        output_fasta=[l]
    return(output_fasta)

def variant_renumber(com_seperated_variantlist,len_whole_seq,strand_correction=True):
    '''takes a comma seperated string of variants and the length of the whole sequence
    renumbers the variants to adjust for the minus strand
    Can be used for 2 types of renumbering: minus strand correction or addition of existing sequence
    in the case of the latter, the "len_whole_seq" actually refers to the len of the existing sequence fragment
    '''
    if com_seperated_variantlist:
        if ',' in com_seperated_variantlist:
            com_seperated_variantlist=com_seperated_variantlist.split(',')
            new_var_pos=[]
            for variant in com_seperated_variantlist:
                if strand_correction:
                    new_var_pos.append(len_whole_seq-int(variant))
                else:
                    new_var_pos.append(len_whole_seq+int(variant))
            return(','.join(new_var_pos))
        else:
            if strand_correction:
                return(str(len_whole_seq-int(com_seperated_variantlist)))
            else:
                return(str(len_whole_seq+int(com_seperated_variantlist))
    return(com_seperated_variantlist)

def import_exons(exoncsv):
    '''read in exons and process the minus strand sequences
    '''
    df=pd.read_csv(exoncsv,sep='\t')
    reverse_strand=df[df['strand']=='-']
    reverse_strand['sequence']=reverse_strand['sequence'].apply(make_complement)
    reverse_strand['pos_het']=reverse_strand.apply(lambda x: variant_renumber(x['pos_het'],x['sequence'].str.len()))
    reverse_strand['pos_hom']=reverse_strand.apply(lambda x: variant_renumber(x['pos_hom'],x['sequence'].str.len()))
    # reverse_strand['pos_het']=reverse_strand['sequence'].str.len()-reverse_strand['pos_het']
    # reverse_strand['pos_hom']=reverse_strand['sequence'].str.len()-reverse_strand['pos_hom']
    return(pd.concat([df[df['strand']!='-'],reverse_strand],axis=0, ignore_index=True))

def make_complement(sequence):
    '''helper function to combine_exons
    makes complimenary mRNA in case of antisense strand'''
    seq=Seq(sequence)
    return(str(seq.reverse_complement()))

def read_gff3_into_frame(gff):
    '''helper function to combine_exons
    input: junction file (gff), and T/F depending on whether you are actually putting in a 0-based (bed) file. True for bed, False for gff.
    this function assumes that the bed file chr is in the foormat "chr2" while gff is in format "2"
    '''
    df=read_csv(gff,sep='\t')
    if gff[-3:]=='bed':
        df=df.drop(df[[3,4,6,7,8]],axis=1) #maybe these shouldn't be hardcoded...
        df.columns=['chromosome','start','end','sense','transcript']
        df['start']=df['start'].astype(int)+1
        df['transcript']=df['transcript'].str.split('transcript:').apply(lambda x: x[1])
        df['transcript'],df['rank']=df['transcript'].str.split('rank=')
        df['transcript']=df['transcript'].str.split(';',1).apply(lambda x: x[0])
        df['rank']=df['rank'].str.split(';',1).apply(lambda x: x[0])
    elif gff[-4:]=='gff3':
        df=df.drop(df[[1,2,5,7]],axis=1)
        df.columns=['chromosome','start','end','sense','transcript']
        df['transcript']=df['transcript'].str.split('transcript_id=').apply(lambda x: x[1])
        df['transcript'],df['rank']=df['transcript'].str.split('exon_number=')
        df['transcript']=df['transcript'].str.split(';',1).apply(lambda x: x[0])
        df['rank']=df['rank'].str.split(';',1).apply(lambda x: x[0])
    else:
        raise Exception("Wrong junction file type input - must be gff3 or bed file")
    df=df[df['chromosome']!='chrM']
    # gff_dict={'chromosome':[],'start':[],'end':[],'sense':[],'transcript':[],'rank':[]}
    # with open(gff) as g:
    #     for line in g:
    #         if 'transcript' in line:
    #             info=re.split('\t',line.strip())
    #             if isBed:
    #                 info_temp=[info[0],info[3],info[4],info[1],info[2]]
    #                 info=info_temp+info[4:]
    #                 chrom=info[0]
    #                 gff_dict['start'].append(int(info[3])+1)
    #             else:
    #                 #chrom='chr'+info[0]
    #                 chrom=info[0] #for gencode
    #                 gff_dict['start'].append(int(info[3]))
    #             gff_dict['sense'].append(info[6])
    #             gff_dict['chromosome'].append(chrom)
    #             gff_dict['end'].append(int(info[4]))
    #             annos=info[-1]
    #             if 'transcript_id' in line: #gencode ref
    #                 tns=annos.split('transcript_id=')[1]
    #                 rk=annos.split('exon_number=')[1]
    #             else: #for ensembl ref
    #                 tns=annos.split('transcript:')[1] 
    #                 rk=annos.split('rank=')[1]
    #             gff_dict['transcript'].append(tns.split(';')[0])
    #             gff_dict['rank'].append(int(rk.split(';')[0]))
    return(df)

def child_initialize(_exon_frame, _junction_frame):
     global exon_frame, junction_frame, terms
     exon_frame = _exon_frame
     junction_frame = _junction_frame

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='track variants')
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--cpu',help='CPU number',required=True)
    parser.add_argument('--exons',help='Phased and variant replaced exons, tab seperated from phase_exons.py',required=True)
    parser.add_argument('--jun', help='file transcript junctions (GFF or BED accepted)', required=True)
    parser.add_argument('--out', help='Output fasta', required=True)
    # parser.add_argument('--report', help='Output file report', required=True) #in favor of reporting the variants in the header of the fasta file
    args=vars(parser.parse_args())
    print('Reading in haplotype-specific exons...')
    # exon_frame=read_exons_into_frame(args['exons'])
    exon_frame=import_exons(args['exons'])
    print('Reading in transcript annotations...')
    junction_frame=read_gff3_into_frame(args['jun']) #can put bed file in
    transcripts=junction_frame.transcript.unique()
    print('building transcripts')
    pool=mp.Pool(processes=int(args['cpu']), initializer=child_initialize, initargs=(exon_frame,junction_frame))
    results=[pool.apply_async(worker_process,args=(transcript,)) for transcript in transcripts] #each transcript processed with seperate cpu
    print_results(results,args['out'])
