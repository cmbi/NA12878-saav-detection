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
                    if id in het:
                        fh.writelines(idline+line)
                    else:
                        f.writelines(idline+line)
    return 'Fastas written'
                    
def id_het(customfa):
    '''sequences with heterozygous variants will have 2 sequences with the same id, save these identifiers'''
    het=Set()
    with open(customfa) as handle:
        recorded=Set()
        for line in handle:
            if line.startswith('>'):
                id=re.split('\|',line,1)[0]
                if id in recorded:
                    het.add(id)
                recorded.add(id)
    return het










