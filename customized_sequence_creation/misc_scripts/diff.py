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

def insilico_digest_diff(cpdtref,cpdtcustom):
    '''this script will find peptides that vary between two cpdt files
    prints cpdt files: one for all different peptides and one with only snv peptides
    '''
    ref=read_cpdt(cpdtref)
    custom=read_cpdt(cpdtcustom)
    newall={}
    newsnv={}
    for id,peplist in ref.iteritems():
        idalt=id+' h1'
        if id in custom or idalt in custom:
            if id in custom:
                key=id
            else:
                key=idalt
            c_peplist=custom[key]
            dif=c_peplist.difference(peplist)
            for pep in dif:
                if determine_snv(pep,peplist):
                    if key in newsnv:
                        newsnv[key].add(pep)
                    else:
                        newsnv[key]=Set()
                        newsnv[key].add(pep)
                if key in newall:
                    newall[key].add(pep)
                else:
                    newall[key]=Set()
                    newall[key].add(pep)
    return newall,newsnv

def determine_snv(peptide,list):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: boolean
    '''
    mismatch=0
    for pep in list:
        if len(pep)==len(peptide):
            for idx,aa in enumerate(pep):
                if aa!=peptide[idx]:
                    mismatch+=1
            if mismatch==1:
                return True
    return False

def read_cpdt(cpdt):
    '''{id:Set(peptides)}'''
    cpdt_pep={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=key.split('|')[0]
                if key in cpdt_pep: #2 haplotypes means 2x the same protein
                    key=key+' h1'
                cpdt_pep[key]=Set()
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(':')[0] #remove the probability
                cpdt_pep[key].add(lp)
    return cpdt_pep

def write_cpdt(dict,outfile):
    '''
    input format: {id:Set(peptides)}
    output: new cpdt file
    '''
    f=open(outfile,'w')
    for id,peplist in dict.iteritems():
        f.writelines(id+'\n')
        for pep in peplist:
            f.writelines('\t'+'PEPTIDE '+pep)
    return 'new CPDT file written to '+outfile









