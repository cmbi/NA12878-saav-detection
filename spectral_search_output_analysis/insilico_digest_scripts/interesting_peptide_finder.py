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
        idho=pid+'_h1'
        idhz=pid+'_h0'
        keys=[]
        if idho in custom or idhz in custom:
            keys=[idhz,idho]
        elif pid in custom:
            keys=[pid]
        if len(keys)>0:
            for key in keys:
                if key in custom:
                    pepdict=custom[key]
                    c_peplist=Set(pepdict.keys())
                    r_peplist=Set(peplist.keys())
                    dif=c_peplist.difference(r_peplist) #save all the peptides that are in the custom but not the reference for this particular protein
                    for pep in dif:
                        if len(pep)>2 and isQualified(pep,ref): #need to check if unique peptide
                            prob=pepdict[pep] #fetch probability
                            if determine_snv(pep,r_peplist):
                                if key in newsnv:
                                    newsnv[key][pep]=prob
                                else:
                                    newsnv[key]={pep:prob}
                            if key in newall:
                                newall[key][pep]=prob
                            else:
                                newall[key]={pep:prob}
    return newall,newsnv

def isQualified(peptide,reference_pepdict):
    for pid,peplist in reference_pepdict.iteritems():
        if peptide in peplist:
            return False
    return True

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
    '''{id:{peptide:probability}'''
    cpdt_pep={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=key.split('|')
                if 'ENSP' in key[0]:
                    key=key[1]
                else:
                    key=key[0]
                cpdt_pep[key]={}
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(': ')
                cpdt_pep[key][lp[0]]=lp[1].strip()
    return cpdt_pep

def write_cpdt(d,outfile):
    '''
    input format: {id:{peptide:probability}}
    output: new cpdt file
    '''
    f=open(outfile,'w')
    for id,peplist in d.iteritems():
        f.writelines('>'+id+'\n')
        for pep,prob in peplist.iteritems():
            f.writelines('\t'+'PEPTIDE '+pep+': '+prob+'\n')
    return 'new CPDT file written to '+outfile

alldiffpeps,snvdiffpeps=insilico_digest_diff(sys.argv[1],sys.argv[2])
print write_cpdt(alldiffpeps,sys.argv[3])
print write_cpdt(snvdiffpeps,sys.argv[4])
