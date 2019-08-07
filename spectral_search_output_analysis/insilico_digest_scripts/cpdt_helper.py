#!/usr/bin/env python2

import re,sys
from sets import Set
#import pandas as pd
#from numpy import median
from __builtin__ import True
from collections import Counter


'''

functions to manipulate/analyze cpdt output

add_labels_to_cpdt_output - put back the headers that are lost upon running cpdt
info_extractor - get some general information from the theoretical digest
discriminant_peptide_finder - looks for discriminant peptides


'''

   
        
def add_labels_to_cpdt_output(cpdtout,originialfa,newcpdt):
    '''
    CP-DT output removes the headers of the fasta. this script puts them back in
    the order is the same
    
    input: the original output of cpdt, the original fasta used to create the cpdt file and the output file name
    output: new cpdt file
    '''
    headers=[]
    body=[]
    with open(originialfa) as ofa:
        for line in ofa:
            if line.startswith('>'):
                headers.append(line)
    with open(cpdtout) as cpdt:
        entry=''
        for line in cpdt:
            if line.strip()!='':
                entry+=line
            elif entry!='':
                body.append(entry)
                entry=''
    if len(headers)==len(body):
        with open(newcpdt,'w') as newcp:
            for idx,fah in enumerate(headers):
                newcp.writelines(fah+body[idx])
    else:
        return "fail, unequal amount of headers vs digests"
    return 'CPDT headers added and written to '+newcpdt

def cpdt_reader(cpdtfile):
    '''
    this function will read a cpdt file and save it into a data format that can be used to compare peptides
    
    this script assumes that the headers have already been put in
    
    {fasta_header:[(pep1,probability),(pep2,probability)]}
    
    input: cpdt file
    output: dictionary
    '''
    cpdt_dict={}
    with open(cpdtfile) as handle:
        header=''
        for line in handle:
            if line.startswith('>'):
                header=line
            elif line.strip().startswith('PEPTIDE'):
                info=re.split('\s|:',line.strip())
                peptide,prob=info[1], info[-1]
                if float(prob)>0.15: #filter out the very low probability peptides, this number is arbitrary
                    cpdt_dict[header].append((peptide,prob))
            elif line.strip()!='':
                header+=line
                cpdt_dict[header]=[]
    return cpdt_dict

def seq_reader(seqfile):
    superseq=''
    with open(seqfile) as ref:
        add=False
        for line in ref:
            if line.startswith('>'):
                add=True
            elif add:
                seq=line.strip()+'|'
                superseq+=seq
                add=False
    return superseq

def info_extractor(dictionary_peptides):
    '''
    gets general information from a dictionary of theoretical peptide digest
    '''
    all_peps=Counter()
    for fa,peplist in dictionary_peptides.iteritems():
        for pep in peplist:
            all_peps[pep[0]]+=1
    print str(len(set(all_peps)))+' unique peptides in this theoretical digest out of '+str(sum(all_peps.values()))+' total peptides.'
    print 'each peptide appears on average '+ str(sum(all_peps.values())/len(all_peps.values())) + ' times and median '+str(median(all_peps.values()))+' times.'
    return 0
        

def discriminant_peptide_finder(cpdtfile,ref_cpdt,output):
    '''
    this function goes through the de novo peptide list and prints out the peptides that can discriminate de novo sequences (or seqs of choice)
    (do i also want to print out the statistics?)
    first reads all the full protein sequences off the ref cpdt file into a Set 
    then looks for each peptide from the de novo set inside the list of sequences and print out only the peptides that do not match existing seqs
    the alternative is to check on the peptide level but I assume it shouldn't matter (the method i am implementing is more strict)
    
    this script assumes that the headers were already put in
    
    input: de novo dictionary made from cpdt_reader function, cpdt ref file (doesn't necessarily have to be cpdt, since fragments are not used)
    output: de novo dictionary printed to a .cpdt file
    '''
    #read in cpdt file
    de_novo_dict=cpdt_reader(cpdtfile)
    #make a long string with all the full length proteins in ref
    ens_seqs=seq_reader(ref_cpdt)
    #check which peptides cannot be matched to those proteins
    de_novo_new={}
    not_matched=Counter()
    matched=Counter()
    matched_prob=Counter()
    for fa,peplist in de_novo_dict.iteritems():
        for tuple in peplist:
            if tuple[0] not in ens_seqs:
                matched[fa]+=1
                matched_prob[tuple[1]]+=1
                if fa not in de_novo_new:
                    de_novo_new[fa]=[tuple]
                else:
                    de_novo_new[fa].append(tuple)
            else:
                not_matched[fa]+=1
    #get information
    print str(len(de_novo_new.keys()))+' proteins have discriminative peptides' #how many proteins have discriminative peptides?
    print 'there were a total of '+str(sum(matched.values()))+ ' discriminative peptides' #how many discriminative peptides were there
    print str(len(list(Counter(el for el in matched.elements() if matched[el]>1)))) #how many proteins had more than 1 discriminative peptide
    #what was the average ratio of discriminative peptides to non-discriminative
    problist=map(float, list(matched_prob.elements()))
    print str(median(problist)) +' is the median probability of seeing a discriminative peptide' #what is the median probability of seeing a discriminative peptide
    #how many of the discriminative peptides have a higher than threshold probability of being seen
    #now print dictionary to cpdt file
    return dict_to_cpdt(de_novo_new,output)

def dict_to_cpdt(dictionary,output):
    if len(dictionary.keys())>0:
        with open(output,'w') as f:
            for header,rest in dictionary.iteritems():
                f.writelines(header+'\n')
                for tup in rest:
                    f.writelines(' PEPTIDE '+tup[0]+': '+tup[1]+'\n')
                f.writelines('\n')            
    return 'peptides written to '+output
            
            
#discriminant_peptide_finder(flair_dict, sys.argv[2], 'discriminative_peptides_ref.cpdt')
#info_extractor(flair_dict)   
add_labels_to_cpdt_output(sys.argv[1],sys.argv[2],sys.argv[3])
            
            
            
            
            
            
