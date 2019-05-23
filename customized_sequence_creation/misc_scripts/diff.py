#!/usr/bin/env python2

import re, collections, sys
from sets import Set

'''

After creating the customized sequences and putting them through ORF prediction
Identify proteins that have differences between the reference and the customized sequences

For the gencode: although a file with predicted protein sequences is provided by gencode and this will likely be the same as predicting with ANGEL, 
it is most ideal to compare ANGEL predicted sequences to other ANGEL predicted sequences to rule out the possibility that sequence prediction played 
a part in the differences between the sequences rather than the variants on the genome level. For this reason it is recommended to use ANGEL predicted
sequence files only

'''

def print_diff_lists(reference,custom,outfa,outfahet):
    '''show differences between the reference and the custom predicted protein sequences
    
    categories: those that reference matches one of the haplotypes, reference matches none of the haplotypes, reference does not match single sequence (het vars absent)
    assumes that both files originate from ANGEL, but the gencode and flair sequences have a different identifier format (hence isONT)
    assumes that the * has been removed from the ANGEL peptide sequence predictions
    
    input: reference peptide sequneces (ANGEL output), customized peptide sequences (ANGEL output), output fasta name, output fasta name (heterozygous)
    output: fastas containing the identifier (from the customized sequences), with reference sequence listed first and then customized protein sequence
    '''
    reference_seqs={} #store id:sequence
    with open(reference) as ref:
        for line in ref:
            if line.startswith('>'):
                tid=re.split('\|',line,1)[0]
            elif line.strip!='':
                reference_seqs[tid]=line.strip()
    f=open(outfa,'w')
    fh=open(outfahet,'w')
    het=id_het(custom)
    with open(custom) as cu:
        for line in cu:
            if line.startswith('>'):
                id=re.split('\|',line,1)[0]
                idline=line
            elif line.strip!='':
                if id in reference_seqs: #theoretically should always be true, but ANGEL orf predictor predicted some sequences in the custom set that it didn't in the reference set. This issue is under investigation
                    if reference_seqs[id]!= line.strip():
                        if id in het:
                            fh.writelines(idline+reference_seqs[id]+'\n'+line)
                        else:
                            f.writelines(idline+reference_seqs[id]+'\n'+line)
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

print print_diff_lists(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])








