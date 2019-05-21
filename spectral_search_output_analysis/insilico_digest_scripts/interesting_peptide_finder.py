#!/usr/bin/env python2

import re, collections, sys
from sets import Set

def insilico_digest_diff(cpdtref,cpdtcustom):
    '''this script will find peptides that vary between two cpdt files
    prints cpdt files: one for all different peptides and one with only snv peptides
    '''
    ref=read_cpdt(cpdtref)
    custom=read_cpdt(cpdtcustom)
    newall={}
    newsnv={}
    for pid,peplist in ref.iteritems():
        idalt=pid+' h1'
        keys=[]
        if idalt in custom:
            keys=[pid,idalt]
        elif pid in custom:
            keys=[pid]
        if len(keys)>0:
            for key in keys:
                c_peplist=custom[key]
                if len(keys)==2 and key==pid:
                    key=pid+' h0'
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

def determine_snv(peptide,plist):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: boolean
    '''
    mismatch=0
    for pep in plist:
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

def write_cpdt(d,outfile):
    '''
    input format: {id:Set(peptides)}
    output: new cpdt file
    '''
    f=open(outfile,'w')
    for id,peplist in d.iteritems():
        f.writelines('>'+id+'\n')
        for pep in peplist:
            f.writelines('\t'+'PEPTIDE '+pep+'\n')
    return 'new CPDT file written to '+outfile

alldiffpeps,snvdiffpeps=insilico_digest_diff(sys.argv[1],sys.argv[2])
print write_cpdt(alldiffpeps,sys.argv[3])
print write_cpdt(snvdiffpeps,sys.argv[4])
