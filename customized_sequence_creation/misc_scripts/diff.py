#!/usr/bin/env python2

import re, collections, sys
from sets import Set
from lib2to3.tests.data import different_encoding

'''

After creating the customized sequences and putting them through ORF prediction
Identify proteins that have differences between the reference and the customized sequences
Quantify interesting outcomes on the protein level due to variations on the genome level such as intron insertions, exon skipping, reading frame shifts, alternative start site, etc

For the gencode: although a file with predicted protein sequences is provided by gencode and this will likely be the same as predicting with ANGEL, 
it is most ideal to compare ANGEL predicted sequences to other ANGEL predicted sequences to rule out the possibility that sequence prediction played 
a part in the differences between the sequences rather than the variants on the genome level. For this reason it is recommended to use ANGEL predicted
sequence files only

'''

def print_diff_lists(reference,custom,outfa,outfahet):
    '''show differences between the reference and the custom predicted protein sequences
    categories: those that reference matches one of the haplotypes, reference matches none of the haplotypes, reference does not match single sequence (het vars absent)'''
    altstart=0
    altend=0
    frameshift=0
    SNP=0
    different=open('different.fa','w')
    undef=open('uncat.fa', 'w')
    reference_seqs={} #store id:sequence
    with open(reference) as ref:
        for line in ref:
            if line.startswith('>'):
                tid=re.split('\|',line,1)[0]
            elif line.strip!='':
                reference_seqs[tid]=line.strip()
    f=open(outfa,'w')
    fh=open(outfahet,'w')
    with open(custom) as cu:
        for line in cu:
            if line.startswith('>'):
                id=re.split('\|',line,1)[0]
                het=False
                if '_h' in id:
                    id=id.split('_h')[0]
                    het=True
                idline=line
            elif line.strip!='':
                if id in reference_seqs:
                    if reference_seqs[id]!= line.strip():
                        category=identify_interesting_events(reference_seqs[id], line.strip())
                        if category==''
                        if het:
                            fh.writelines(idline+reference_seqs[id]+'\n'+line)
                        else:
                            f.writelines(idline+reference_seqs[id]+'\n'+line)
    return 'Fastas written'

def identify_interesting_events(reference,custom):
    ''' input reference sequence and custom sequence and classify the event that occurred
    possible categories:
        - alternative start
        - alternative end
        - indel
        - frame shift
        - SNP
        - other different
        - undefined
    use align when necessary?
    '''
    if len(reference)==len(custom):
        consecutive=0
        streak=[]
        mismatch=0
        for idx,aa in enumerate(reference):
            if aa!=custom[idx]:
                consecutive+=1
                mismatch+=1
            else:
                if consecutive>1:
                    streak.append(consecutive)
                consecutive=0
        if consecutive>1:
            streak.append(consecutive)
        if consecutive>15:
            return "frameshift"
        if len(streak)>0 and sum(streak)>20:
            return 'different'
        elif len(streak)==0 and mismatch>0:
            return 'SNP'
    elif len(reference)>len(custom): #custom shorter
        if custom in reference:
            rest=reference.split(custom)
            if rest[0]=='':
                return "alternative end"
            elif rest[1]=='':
                return "alternative start"
        elif len(custom)>15: #this is the tricky part, conservative tag based search
            return tag_based_search(custom, reference)
    elif len(custom)>len(reference): #reference shorter
        if reference in custom:
            rest=custom.split(reference)
            if rest[0]=='':
                return "alternative end"
            elif rest[1]=='':
                return "alternative start"
        elif len(reference)>15:
            return tag_based_search(reference, custom)
    return event

def tag_based_search(shorter,longer):
    '''in the case where imperfect match between 2 varied length strings, use tag based search'''
    if shorter[:6] and shorter[-6:] in longer: #in case the shorter is completely inside the longer
        head=re.split(shorter[:6],longer,1)
        tail=head[1].split(shorter[-6:])
        if len(tail)>2: #if the end pattern 
            tail=[''.join(tail[:-1]),tail[-1]]
        if len(head[0])<3 and len(tail[1])<3: #something in the middle that is different
            return "indel"
        elif len(head[0])>3 and abs(head[1]-len(shorter))<3:
            return "alternative start"
        elif len(tail[-1])>3 and abs(len(head[0])+len(tail[0])-len(shorter))<3:
            return "alternative end"
    return "undefined" #doesn't deal with SNPs in the tags

def het_diff_info(difffile):
    '''check how many sequences had only one haplotype that was different from the reference rather than both
    input: the het file generated from print_diff_lists function'''
    onehap=Set()
    twohap=Set()
    with open(difffile) as handle:
        for line in handle:
            if line.startswith('>'):
                tid=re.split('_h',line,1)[0]
                if tid in onehap:
                    onehap.remove(tid)
                    twohap.add(tid)
                else:
                    onehap.add(tid)
    return str(len(onehap))+' sequences had one haplotype differing from the reference while '+str(len(twohap))+' had both haplotypes differing from the reference'
    
print print_diff_lists(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])








