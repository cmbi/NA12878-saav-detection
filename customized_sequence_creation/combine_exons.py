#!/usr/bin/env python3

import re, collections, sys
import pandas as pd
import numpy as np
#from enum import Enum
from Bio.Seq import Seq
from tqdm import tqdm
import multiprocessing as mp

'''
combine exons: given "phased" exon sequences from phase_exons, put the exons back together (per haplotype) to create full length transcripts
    -read_exons_into_frame: parses the exon file and organizes them into a pandas dataframe
    -read_gff_into_frame: parses the junction files and organizes them into a pandas dataframe
    -make_complement: in case minus strand, make reverse complement
'''
                        
def print_results(results,outfile,routfile):
    '''
    takes the exons created from the phase_exons function and stitches them together according to the gff file to make full transcripts
    
    ''' 
    f=open(outfile,'w')
    r=open(routfile,'w')
    output = [p.get() for p in results] #list of tuplies
    #write transcript to fasta
    for o in output:
        transcript=o[0]
        for l in o[1]:
            f.write(l)
        var=o[2]
        org=o[3]
        if len(o[1])>1:
            r.writelines(transcript+'_h0\t'+','.join(var[0])+'\t'+','.join(org[0])+'\n')
            r.writelines(transcript+'_h1\t'+','.join(var[1])+'\t'+','.join(org[1])+'\n')
        elif len(var[0])>0:
            r.writelines(transcript+'\t'+','.join(var[0])+'\t'+','.join(org[0])+'\n')
    return('done')

def worker_process(transcript):
    sequence_zero=''
    sequence_one=''
    to_add=''
    kill=False
    AS=False
    variants=[[],[]]
    origin=[[],[]]
    to_write=[]
    trans_info=junction_frame[junction_frame["transcript"]==transcript] #get junction information for the transcript
    trans_info=trans_info.sort_values(by="rank") #sort them in the correct order
    for tidx,trow in trans_info.iterrows(): #iterate through each junction and get the sequence
        sense=trow["sense"]
        seq_rows=exon_frame[(exon_frame["chromosome"]==trow["chromosome"]) & (exon_frame["end"]==trow["end"]) & (exon_frame["start"]==trow["start"])]
        if len(seq_rows.index)==1: #if no heterozygous variants in exon
            #gather variants
            v=seq_rows.iloc[0]["position"]
            v=process_variant(v)
            o=seq_rows.iloc[0]["origin"]
            o=process_origin(o)
            if sense=='-':
                to_add=make_complement(seq_rows.iloc[0]["sequence"])
                if v!='NA':
                    v=[len(to_add)-x for x in v]
            else:
                to_add=seq_rows.iloc[0]["sequence"]
            len_existingz,len_existingo=len(sequence_zero),len(sequence_one)
            if v!='NA':
                variants[0].extend([str(x+len_existingz) for x in v])
                variants[1].extend([str(x+len_existingo) for x in v])
                origin[0].extend(o)
                origin[1].extend(o)
            sequence_zero+=to_add
            sequence_one+=to_add
        elif len(seq_rows.index)==2: #if heterozygous variants in exon
            AS=True
            seq_row_hz=seq_rows[seq_rows["haplotype"]==0]
            seq_row_ho=seq_rows[seq_rows["haplotype"]==1]
            vz=seq_row_hz.iloc[0]["position"]
            vo=seq_row_ho.iloc[0]["position"]
            vz=process_variant(vz)
            vo=process_variant(vo)
            oo=seq_row_ho.iloc[0]["origin"]
            oz=seq_row_hz.iloc[0]["origin"]
            oo=process_origin(oo)
            oz=process_origin(oz)
            to_add_z=make_complement(seq_row_hz.iloc[0]["sequence"])
            to_add_o=make_complement(seq_row_ho.iloc[0]["sequence"])
            len_existingz,len_existingo=len(sequence_zero),len(sequence_one)
            if sense=='-':
                sequence_zero+=make_complement(to_add_z)
                sequence_one+=make_complement(to_add_o)
                if vz!='NA':
                    vz=[len(to_add_z)-x for x in vz]
                if vo!='NA':
                    vo=[len(to_add_o)-x for x in vo]
            else:
                sequence_zero+=seq_row_hz.iloc[0]["sequence"]
                sequence_one+=seq_row_ho.iloc[0]["sequence"]
            if vz!='NA' and vo!='NA':
                variants[0].extend([str(x+len_existingz) for x in vz])
                variants[1].extend([str(x+len_existingo) for x in vo])
                origin[0].extend(oz)
                origin[1].extend(oo)
        else:
            print('exon not found')
            print(trans_info)
            kill=True
            sys.exit()
    if not kill:
        if AS:
            # count_allele_specific+=1
            l1='>'+transcript+'_h0'+'\n'+sequence_zero+'\n'
            l2='>'+transcript+'_h1'+'\n'+sequence_one+'\n'
            to_write=[l1,l2]
        else:
            l='>'+transcript+'\n'+sequence_zero+'\n'
            to_write=[l]
    return([transcript,to_write,variants,origin])

def process_variant(v):
    if v!='NA':
        if ',' in v:
            v=v.split(',')
            v=map(int,v)
        else:
            v=[int(v)]
    return(v)

def process_origin(o):
    if o!='NA':
        if ',' in o:
            o=o.split(',')
        else:
            o=[o]
    return(o)

def make_complement(sequence):
    '''helper function to combine_exons
    makes complimenary mRNA in case of antisense strand'''
    seq=Seq(sequence)
    return(str(seq.reverse_complement()))

def read_exons_into_frame(hap_exons):
    '''helper function to combine_exons
    reads exon file and stores into data frame
    this function contains the correction for the 0/1 coordinates in variable "true_start"
    '''                
    exon_dict={'chromosome':[],'start':[],'end':[],'haplotype':[],'position':[],'origin':[],'sequence':[]}
    with open(hap_exons) as he:
        for line in he:
            if line.startswith('>'):
                if 'chrY' not in line and 'chrM' not in line: #why did i block out the M chromosome again?
                    if 'haplotype' in line:
                        info=line.strip().split(' ')
                        hap=info[1].split(':')[1]
                        pos=info[2].split(':')[1]
                        org=info[3].split(':')[1]
                        exon_dict['haplotype'].append(int(hap))
                        exon_dict['position'].append(pos)
                        exon_dict['origin'].append(org)
                        info=re.split(':|-',info[0])
                    elif 'pos' in line:
                        info=line.strip().split(' ')
                        pos=info[1].split(':')[1]
                        org=info[2].split(':')[1]
                        exon_dict['position'].append(pos)
                        exon_dict['origin'].append(org)
                        exon_dict['haplotype'].append(9)
                        info=re.split(':|-',info[0])
                    else:
                        exon_dict['haplotype'].append(9)
                        exon_dict['position'].append('NA')
                        exon_dict['origin'].append('NA')
                        info=re.split(':|-',line.strip())
                    chrom=info[0]
                    exon_dict['chromosome'].append(chrom[1:])
                    true_start=int(info[1])+1
                    exon_dict['start'].append(true_start)
                    exon_dict['end'].append(int(info[2]))
            else:
                exon_dict['sequence'].append(line.strip())
    return(pd.DataFrame(data=exon_dict))

def read_gff3_into_frame(gff):
    '''helper function to combine_exons
    input: junction file (gff), and T/F depending on whether you are actually putting in a 0-based (bed) file. True for bed, False for gff.
    this function assumes that the bed file chr is in the foormat "chr2" while gff is in format "2"
    '''
    if gff[-3:]=='bed':
        bed=True
    else:
        bed=False
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
    return(pd.DataFrame(data=gff_dict))

def child_initialize(_exon_frame, _junction_frame):
     global exon_frame, junction_frame, terms
     exon_frame = _exon_frame
     junction_frame = _junction_frame

if __name__ == '__main__':
    hap_exons=sys.argv[1]
    gff=sys.argv[2]
    outfile=sys.argv[3]
    routfile=sys.argv[4]
    print('Reading in haplotype-specific exons...')
    exon_frame=read_exons_into_frame(hap_exons)
    print('Reading in transcript annotations...')
    junction_frame=read_gff3_into_frame(gff) #can put bed file in
    transcripts=junction_frame.transcript.unique()
    killed=0
    print('building transcripts')
    pool_size=45
    pool=mp.Pool(processes=pool_size, initializer=child_initialize, initargs=(exon_frame,junction_frame))
    results=[pool.apply_async(worker_process,args=(transcript,)) for transcript in transcripts]
    print_results(results,outfile,routfile)
