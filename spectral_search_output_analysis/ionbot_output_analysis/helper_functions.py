#!/usr/bin/env python3

import re, os,sys, itertools
import pandas as pd
from collections import Counter



def find_chrom(prots,chromdict):
    for p in prots:
        p=get_id(p)
        if '_h' in p:
            p=p.split('_h')[0]
        if p in chromdict:
            # return(re.sub('r','r ',chromdict[p]))
            chrom=chromdict[p].split('chr')[1]
            if chrom.isdigit():
                return(int(chrom))
            else:
                return(chrom)
    return("unknown")

def find_strand(prots,stranddict):
    for p in prots:
        p=get_id(p)
        if '_h' in p:
            p=p.split('_h')[0]
        if p in stranddict:
            if stranddict[p]=='-':
                return("reverse")
            elif stranddict[p]=='+':
                return("forward")
            # return(stranddict[p])
    return("unknown")

def get_id(idstring):
    i=idstring
    if '|m.' in i:
        i=i.split('|m.')[0]
    elif 'ENSP' in i:
        i=i.split('|')[1]
    return(i)

def bin_hits_by_source(scanid,ids,oro,ono,ob,isOpenmut):
    '''sort peptide hits by their source dictionary'''
    ref_only=oro
    ont_only=ono
    both=ob
    ont=False
    ref=False
    for i in ids:
        if not isOpenmut:
            if '_E' in i:
                ont=True
            elif 'ENST' in i:
                ref=True
        else:
            if '|m.' in i: #assumes that the proteins from the reference were not generated with ANGEL
                ont=True
            elif 'ENST' in i:
                ref=True
    if ont and ref:
        both.add(scanid)
    elif ont:
        ont_only.add(scanid)
    elif ref:
        ref_only.add(scanid)
    #else:
    #    raise Exception('Unexpected protein found')
    return(ref_only,ont_only,both)


def counter_translator(counterobj):
    ct_prots=Counter()
    for elem,ct in counterobj.items():
        if ct>19:
            ct_prots['20+']+=1
        else:
            ct_prots[ct]+=1
    return(ct_prots)

def fill_cpdt(pep,ids,old_cpdt_pep):
    '''fill the data structure to match predicted (mutated) peptides to observed
    
    how many unique variant peptides are detected, from how many unique proteins?
    how many total instances of correct/incorrect variant peptides are detected?
    build up the dictionary with every iteration of the 
    '''
    cpdt_pep=old_cpdt_pep
    notfound=0
    for isi in ids:
        i=get_id(isi)
        found=False
        possibilities=[i,i+'_h0',i+'_h1']
        for poss in possibilities:
            if poss in cpdt_pep.keys():
                ispeplist=cpdt_pep[poss]#make copy to mutate as you iterate
                for p,ct in cpdt_pep[poss].items():
                    if p==pep:
                        ct+=1
                        ispeplist[p]=ct
                        found=True
                cpdt_pep[poss]=ispeplist
                # break #only counting the first of the 2 haplotypes- my choices are to pick one or double count, better to double count?
            if not found:
                notfound+=1
                ##this recovers a lot of peptides that would otherwise be filtered out
                # if '->' not in mod and found==False:
                #     if pep in full_seqs[i] or 'X' in full_seqs[i]:
                #         ispeplist[pep]=1
                #         cpdt_pep[i]=ispeplist
    if notfound==len(ids):
        missed=1
    else:
        missed=0
    return(cpdt_pep,missed)

def initiate_counter():
    '''generate all possible AA subsititutions and put them in a counter'''
    all_aa=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    all_list=list(itertools.product(all_aa,all_aa))
    all_counter=Counter()
    for l in all_list:
        all_counter[l]=0
    return(all_counter)

# def counter_to_df(allc, observed):
#     '''create the df that will be used in the heatmap figure, use normalization (min-max)'''
#     df=get_normalized_matrix(allc,observed)
#     return(df,dfall)

def determine_snv(peptide,plist):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: boolean
    '''
    for pep in plist:
        if len(pep)==len(peptide):
            original=''
            sub=''
            mismatch=0
            for idx,aa in enumerate(pep):
                if aa!=peptide[idx]:
                    original=aa
                    sub=peptide[idx]
                    situation1= original=='I' and sub=='L'
                    situation2= original=='L' and sub=='I'
                    if not situation1 and not situation2:
                        mismatch+=1
            if mismatch==1:
                return(pep,(original,sub))
    return('','')

def detect_peptides(pep,ids,cpdt_pep,isOpenmut):
    for isi in ids:
        i=get_id(isi)
        found=False
        if isOpenmut:
            possibilities=[i,i+'_h0',i+'_h1']
        else:
            possibilities=[i]
        for poss in possibilities:
            if poss in cpdt_pep.keys():
                for p in cpdt_pep[poss]: #this allows for some flexibility in matches, allows for non-canonical
                    if pep in p or equivalent_check(pep,p):
                        return(poss,p)
    return('','')

def equivalent_check(querypep,pep):
    chunks=re.split('I|L',pep)
    chunks_query=re.split('I|L',querypep)
    if contains(chunks_query,chunks):
        return(True)
    return(False)

def contains(small, big):
    for i in range(1 + len(big) - len(small)):
        if small == big[i:i+len(small)]:
            return(True)
    return(False)

def add_to_observed(pep,ids,cpdtpep,isOpenmut):
    mut_prot,mut_pep=detect_peptides(pep,ids,cpdtpep,isOpenmut)
    #add mutant peptide to observed
    newdict=cpdtpep
    if mut_prot!='':
        # if mut_prot not in newdict:
        #     newdict[mut_prot]=Counter()
        newdict[mut_prot][mut_pep]+=1
    return(newdict)


def get_normalized_matrix(allc,observed):
    df_all = pd.DataFrame(list(allc.values()),index=pd.MultiIndex.from_tuples(allc.keys()),columns=['values'])
    df_observed = pd.DataFrame(list(observed.values()),index=pd.MultiIndex.from_tuples(observed.keys()),columns=['values'])
    df_all['normalized']=(df_all['values']-df_all['values'].min())/(df_all['values'].max()-df_all['values'].min()) # if want to do min-max normalization
    df_observed['normalized']=(df_observed['values']-df_observed['values'].min())/(df_observed['values'].max()-df_observed['values'].min()) # if want to do min-max normalization
    df_combi=df_observed['normalized']-df_all['normalized']
    #df['normalized']=(df['values']-df['values'].mean())/df['values'].std() # if want to do standard normalization
    ser=pd.Series(df_combi)
    df = ser.unstack()#.fillna(0)
    return(df)

def count_muts(full_cpdt_dict):
    allmuts=Counter()
    for prot,mutct in full_cpdt_dict.items(): #no open mutation
        for pepi,pepct in mutct.items():
            if pepct>0:
                allmuts[pepi]+=pepct
    return(allmuts)

def detected_proteins(ids,pco):
    proteins_covered=pco
    for idu in ids:
        proteins_covered[idu]+=1
    return(proteins_covered)

def categorize_mods(list_mods):
    mod_ct=Counter()
    for mod in list_mods:
        mod=str(mod)
        if len(re.findall('[A-Z]->[A-Z]',mod))>0:
            mod_ct['SAAV']+=1
        elif mod=='nan':
            mod_ct["none"]+=1
        else:
            s_mod=re.split('\[[a-z]\]',mod)[0]
            mod_ct[s_mod]+=1
    return(mod_ct)

def gather_counts(peptide_counter):
    length_counter=Counter()
    for pep in peptide_counter:
        length_counter[len(pep)]+=peptide_counter[pep]
    return(length_counter)
