#!/usr/bin/env python2

import re

'''bed12_to_bed: converts a bed12 file to a 1 exon per line bed file, to then be used with getfasta'''

def bed12_to_bed(bed,newbed):
    '''change bed12 to a 1-line-per-exon format that is easy for bedtools getfasta to make into a fasta file
    This is necessary for the flair fasta file since there is no gff3 for it (only psl, which was then converted into bed12 with bedtools)
    for every line of the bed file
        transfer chromosome over directly
        iterate through every start position, add length to make end position
        print to line with transcript name and rank
    '''
    f=open(newbed,'w')
    with open(bed) as b:
        for line in b:
            l=line.strip().split('\t')
            chr=l[0]
            name=l[3]
            strand=l[5]
            starts=l[-1].split(',')[:-1]
            lengths=l[-3].split(',')[:-1]
            ends=[]
            for i,s in enumerate(starts):
                info='transcript:'+name+';rank='+str(i+1)+';'
                end=int(s)+int(lengths[i])
                f.writelines('\t'.join([chr,s,str(end),'.','.','.',strand,'.','.',info])+'\n')


bed12_to_bed(sys.argv[1],sys.argv[2])