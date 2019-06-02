#!/usr/bin/env python2

import re,sys
from sets import Set
import pandas as pd
from numpy import median
from __builtin__ import True
from collections import Counter


'''
this script will contain functions for peptide dictionary manipulation

goal1: REFERENCE: have a complete dictionary of "reference", with all ref annotated proteins combined with orf predictions for non-protein coding transcripts as well as orf prediction for novel flair isoforms
goal2: CUSTOM: using the variant inclusive transcript ORF predictions


one potential issue: if replacing a sequence in the reference, should we use the reference sequence orf to update a (potentially incorrect) predicted orf from the cell line transcripts?

'''


def cpdt_reader(cpdtfile):
    '''
    this function will read a cpdt file and save it into a data format that can be used to compare peptides
    
    this script assumes that the headers have already been put in
    
    {fasta_header:[(pep1,probability),(pep2,probability)]}
    
    input: cpdt file
    output: dictionary
    '''
    cpdt_dict={}
    all_peps=Counter()
    with open(cpdtfile) as handle:
        header=''
        for line in handle:
            if line.startswith('>'):
                header=line
            elif line.strip().startswith('PEPTIDE'):
                info=re.split('\s|:',line.strip())
                peptide,prob=info[1], info[-1]
                all_peps[peptide]+=1
                if float(prob)>0.15: #filter out the very low probability peptides, this number is arbitrary
                    cpdt_dict[header].append((peptide,prob))
            elif line.strip()!='':
                header+=line
                cpdt_dict[header]=[]
    print str(len(set(all_peps)))+' unique peptides in this theoretical digest out of '+str(sum(all_peps.values()))+' total peptides.'
    print 'each peptide appears on average '+ str(sum(all_peps.values())/len(all_peps.values())) + ' times and median '+str(median(all_peps.values()))+' times.'
    return cpdt_dict
        

def discriminant_peptide_finder(transcriptome_cpdt,ref_cpdt,output):
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
    de_novo_dict=cpdt_reader(transcriptome_cpdt)
    info_extractor(de_novo_dict)
    #first make a long string with all the full length proteins in ref
    ens_seqs=''
    with open(ref_cpdt) as ref:
        add=False
        for line in ref:
            if line.startswith('>'):
                add=True
            elif add:
                seq=line.strip()+'|'
                ens_seqs+=seq
                add=False
    #then check which peptides cannot be matched to those proteins
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
            
transcriptome_cpdt=sys.argv[1]
reference_cpdt=sys.argv[2]
discriminant_peptide_finder(transcriptome_cpdt, reference_cpdt, 'discriminative_peptides_ref.cpdt')
            
            
            
            
            
            
