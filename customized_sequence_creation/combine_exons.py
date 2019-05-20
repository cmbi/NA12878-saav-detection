#!/usr/bin/env python2

import re, collections, sys
from sets import Set
import pandas as pd
import numpy as np
#from enum import Enum
from Bio.Seq import Seq
from tqdm import tqdm

'''
combine exons: given "phased" exon sequences from phase_exons, put the exons back together (per haplotype) to create full length transcripts
    -read_exons_into_frame: parses the exon file and organizes them into a pandas dataframe
    -read_gff_into_frame: parses the junction files and organizes them into a pandas dataframe
    -make_complement: in case minus strand, make reverse complement
'''
                        
def combine_exons(hap_exons,gff,isBed,outfile):
    '''
    takes the exons created from the phase_exons function and stitches them together according to the gff file to make full transcripts
    
    ''' 
    f=open(outfile,'w')
    sequence_zero=''
    sequence_one=''
    to_add=''
    count_allele_specific=0
    print 'Reading in haplotype-specific exons...'
    exon_frame=read_exons_into_frame(hap_exons)
    print 'Reading in transcript annotations...'
    junction_frame=read_gff3_into_frame(gff,isBed) #add true (if bed) or false (if gff) here
    transcripts=junction_frame.transcript.unique()
    killed=0
    print 'building transcripts'
    for transcript in tqdm(transcripts): #iterate through unique transcript ids
        kill=False
        AS=False
        trans_info=junction_frame[junction_frame["transcript"]==transcript] #get junction information for the transcript
        trans_info=trans_info.sort_values(by="rank") #sort them in the correct order
        for tidx,trow in trans_info.iterrows(): #iterate through each junction and get the sequence
            sense=trow["sense"]
            seq_rows=exon_frame[(exon_frame["chromosome"]==trow["chromosome"]) & (exon_frame["end"]==trow["end"]) & (exon_frame["start"]==trow["start"])]
            if len(seq_rows.index)==1: #if no heterozygous variants in exon
                if sense=='-':
                    to_add=make_complement(seq_rows.iloc[0]["sequence"])
                else:
                    to_add=seq_rows.iloc[0]["sequence"]
                sequence_zero+=to_add
                sequence_one+=to_add
            elif len(seq_rows.index)==2: #if heterozygous variants in exon
                AS=True
                seq_row_hz=seq_rows[seq_rows["haplotype"]==0]
                seq_row_ho=seq_rows[seq_rows["haplotype"]==1]
                if sense=='-':
                    sequence_zero+=make_complement(seq_row_hz.iloc[0]["sequence"])
                    sequence_one+=make_complement(seq_row_ho.iloc[0]["sequence"])
                else:
                    sequence_zero+=seq_row_hz.iloc[0]["sequence"]
                    sequence_one+=seq_row_ho.iloc[0]["sequence"]
            else:
                print 'exon not found'
                print trans_info
                kill=True
                killed+=1
        #write transcript to fasta
        if not kill:
            if AS:
                count_allele_specific+=1
                f.writelines('>'+transcript+' haplotype:0'+'\n'+sequence_zero+'\n')
                f.writelines('>'+transcript+' haplotype:1'+'\n'+sequence_one+'\n')
            else:
                f.writelines('>'+transcript+'\n'+sequence_zero+'\n')
        sequence_zero=''
        sequence_one=''
    return 'number of allele specific transcripts: '+str(count_allele_specific)+ '. number exons not found= '+str(killed)

def make_complement(sequence):
    '''helper function to combine_exons
    makes complimenary mRNA in case of antisense strand'''
    seq=Seq(sequence)
    return str(seq.reverse_complement())

def read_exons_into_frame(hap_exons):
    '''helper function to combine_exons
    reads exon file and stores into data frame
    this function contains the correction for the 0/1 coordinates in variable "true_start"
    '''                
    exon_dict={'chromosome':[],'start':[],'end':[],'haplotype':[],'sequence':[]}
    with open(hap_exons) as he:
        for line in he:
            if line.startswith('>'):
                if 'chrY' not in line and 'chrM' not in line:
                    if 'haplotype' in line:
                        info=line.strip().split(' ')
                        hap=info[1].split(':')[1]
                        exon_dict['haplotype'].append(int(hap))
                        info=re.split(':|-',info[0])
                    else:
                        exon_dict['haplotype'].append(9)
                        info=re.split(':|-',line.strip())
                    chrom=info[0]
                    exon_dict['chromosome'].append(chrom[1:])
                    true_start=int(info[1])+1
                    exon_dict['start'].append(true_start)
                    exon_dict['end'].append(int(info[2]))
            else:
                exon_dict['sequence'].append(line.strip())
    return pd.DataFrame(data=exon_dict)

def read_gff3_into_frame(gff,bed):
    '''helper function to combine_exons
    input: junction file (gff), and T/F depending on whether you are actually putting in a 0-based (bed) file. True for bed, False for gff.
    this function assumes that the bed file chr is in the foormat "chr2" while gff is in format "2"
    '''
    gff_dict={'chromosome':[],'start':[],'end':[],'sense':[],'transcript':[],'rank':[]}
    with open(gff) as g:
        for line in g:
            if 'transcript' in line:
                info=re.split('\t',line.strip())
                if bed:
                    info_temp=[info[0],info[3],info[4],info[1],info[2]]
                    info=info_temp+info[4:]
                    chrom=info[0]
                    gff_dict['start'].append(int(info[3])+1)
                else:
                    #chrom='chr'+info[0]
                    chrom=info[0] #for gencode
                    gff_dict['start'].append(int(info[3]))
                gff_dict['sense'].append(info[6])
                gff_dict['chromosome'].append(chrom)
                gff_dict['end'].append(int(info[4]))
                annos=info[-1]
                if 'transcript_id' in line: #gencode ref
                    tns=annos.split('transcript_id=')[1]
                    rk=annos.split('exon_number=')[1]
                else: #for ensembl ref
                    tns=annos.split('transcript:')[1] 
                    rk=annos.split('rank=')[1]
                gff_dict['transcript'].append(tns.split(';')[0])
                gff_dict['rank'].append(int(rk.split(';')[0]))
    return pd.DataFrame(data=gff_dict)
    

print combine_exons(sys.argv[1], sys.argv[2],sys.argv[3])
