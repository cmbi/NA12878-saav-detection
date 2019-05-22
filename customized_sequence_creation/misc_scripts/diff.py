#!/usr/bin/env python2

import re, collections, sys
from sets import Set

'''Identify proteins that have differences between the reference and the customized sequences'''

def print_diff_lists(reference,custom,outfa,outfahet):
    '''show differences between the reference and the custom predicted protein sequences
    categories: those that reference matches one of the haplotypes, reference matches none of the haplotypes, reference does not match single sequence (het vars absent)'''
    reference_seqs={} #store id:sequence
    with open(reference) as ref:
        for line in ref:
            if line.startswith('>'):
                id=re.split('\|',line,1)[0]
            elif line.strip!='':
                reference_seqs[id]=line.strip()
    f=open(outfa,'w')
    fh=open(outfahet,'w')
    het=id_het(custom)
    with open(custom) as cu:
        for line in cu:
            if line.startswith('>'):
                id=re.split('\|',line,1)[0]
                idline=line
            elif line.strip!='':
                if reference_seqs[id]!= line.strip():
                    idh=id+' h1'
                    if id in het or idh in het:
                        fh.writelines(idline+line)
                    else:
                        f.writelines(idline+line)
    return 'Fastas written'
                    
def id_het(customfa):
    het=Set()
    with open(customfa) as handle:
        previd=''
        for line in handle:
            if line.startswith('>'):
                id=re.split('\|',line,1)[0]
                if id==previd:
                    het.add(id)
                previd=id
    return het










