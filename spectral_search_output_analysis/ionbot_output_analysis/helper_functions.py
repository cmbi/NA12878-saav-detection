#!/usr/bin/env python3

import re, os,sys, itertools
import pandas as pd
from collections import Counter
import file_import



def find_chrom(prots,chromdict):
    for p in prots:
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
        if '_h' in p:
            p=p.split('_h')[0]
        if p in stranddict:
            if stranddict[p]=='-':
                return("reverse")
            elif stranddict[p]=='+':
                return("forward")
            # return(stranddict[p])
    return("unknown")

def get_id(idstring,base=False):
    if '|m.' in idstring:
        outstring=idstring.split('|m.')[0]
    elif 'ENSP' in idstring:
        if '|' not in idstring:
            return('none')
        tid=idstring.split('|')[1]
        if 'Random' in idstring:
            prefix=idstring.split('_')[0]
            outstring=prefix+'_'+tid
        else:
            outstring=tid
    else:
        outstring=idstring.strip()
    if base and '_h' in outstring:
        return(outstring.split('_h')[0])
    return(outstring)

def bin_hits_by_source(protstring):
    '''sort peptide hits by their source dictionary'''
    sources=set()
    if 'Random' in protstring or 'random' in protstring:
        return('decoy')
    if '||' in protstring:
        ids=protstring.split('||')
    else:
        ids=[protstring]
    for i in ids:
        if '_ENSG' in i:
            sources.add('ont')
        else:
            sources.add('ref')
    if 'ont' in sources and 'ref' in sources:
        return('both')
    return(list(sources)[0])


def counter_translator(counterobj):
    ct_prots=Counter()
    for elem,ct in counterobj.items():
        if ct>19:
            ct_prots['20+']+=1
        else:
            ct_prots[ct]+=1
    return(ct_prots)

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

def determine_snv(peptide,counterpart):
    ''' extract the single amino acid substitution from the 
    input: a peptide and a list of peptides
    output: tuple
    '''
    assert len(peptide)==len(counterpart), f"fatal error- counterpart and variant peptide length don't match"
    for a, b in zip(peptide, counterpart):
        if a != b and a!='*' and b!='*':
            return(b,a)
    return(None)

def detect_peptides(pep,ids,cpdt_pep,isOpenmut,debug=False,include_extra=False):
    count=0
    for isi in ids:
        found=False
        if isOpenmut:
            possibilities=[isi,isi+'_h0',isi+'_h1']
        else:
            possibilities=[isi]
        for poss in possibilities:
            if poss in cpdt_pep.keys():
                for p in cpdt_pep[poss]: #this allows for some flexibility in matches, allows for non-canonical
                    if pep in p or equivalent_check(pep,p):
                        if debug:
                            # count+=1
                            print(pep,cpdt_pep[poss])
                        if include_extra:
                            return(True,p,poss)
                        return(True)
    if include_extra:
        return(False,'','')
    return(False)

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

def abbreviate_peps(counter_varpep):
    '''
    change the counter to only have one peptide per SAAV to simplify analysis later on
    '''
    set_varpep=set(counter_varpep)
    vp=sorted(set_varpep,key=len) #short to long
    accounted_for=set()
    abbreviated=Counter()
    abbr_dict={}
    for p in vp:
        if p not in accounted_for:
            rest=set_varpep.difference(accounted_for)
            rest.discard(p)
            added=False
            for r in rest:
                if len(p)>len(r) and contains(r,p):
                    print('something wrong')
                elif len(r)>len(p) and contains(p,r):
                    if p in abbreviated:
                        abbreviated[p]+=counter_varpep[r]
                        abbr_dict[p].append(r)
                    else:
                        abbreviated[p]=counter_varpep[r]+counter_varpep[p]
                        abbr_dict[p]=[p,r]
                    accounted_for.add(r)
                    accounted_for.add(p)
                    added=True
            if not added:
                abbreviated[p]=counter_varpep[p]
                abbr_dict[p]=[p]
                accounted_for.add(p)
                # for a in abbreviated:
                #     if len(a)>len(r) and contains(r,a):
                #         adtn_ct=counter_varpep[a]
    return(abbreviated,abbr_dict)

def unpack_grouped_variants(set_shortpeps,reference_dict):
    ''' return set of all possible peptides with a particular variant given a set of the shortest peptides for a set of variants
    '''
    concat=set()
    for pep in set_shortpeps:
        if pep in reference_dict:
            concat=concat.union(set(reference_dict[pep]))
    return(concat)

def concat_dicts(x,y):
    c={}
    all_keys=set(x.keys()).union(set(y.keys()))
    for a in all_keys:
        if a in x and a in y:
            c[a]=set(x[a]).union(set(y[a]))
        elif a in x:
            c[a]=set(x[a])
        else:
            c[a]=set(y[a])
    return(c)

def add_to_observed(pep,ids,cpdtpep,isOpenmut,debug=False):
    mut,mut_pep,mut_prot=detect_peptides(pep,ids,cpdtpep,isOpenmut,include_extra=True)
    #add mutant peptide to observed
    newdict=cpdtpep
    if mut:
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

def categorize_mods(mod):
    if len(re.findall('[A-Z]->[A-Z]',mod))>0:
        return('SAAV')
    elif mod=='nan':
        return("none")
    s_mod=re.split('\[[a-z]\]',mod)[0]
    return(s_mod)

def take_first_protein(protstring):
    if '||' in protstring:
        return(get_id(protstring.split('||')[0],base=True))
    return(get_id(protstring,base=True))

def remove_empty(variant_pep_dict):
    '''remove empty entries from the variant peptide dictionary'''
    newdict={}
    for prot,peplist in variant_pep_dict.items():
        for pep,ct in peplist.items():
            if ct>0:
                if prot not in newdict:
                    newdict[prot]={}
                newdict[prot][pep]=ct
    return(newdict)

def longest(s):
    return(max(s,key=len))

def get_all_observed(vf,vc,theoretical):
    all_peptides=file_import.il_sensitive_read_csv(theoretical,names=['protein','peptide','start'],to_replace=['peptide'],variant=False)
    return(vf.merge(all_peptides, on='peptide'),vc.merge(all_peptides, on='peptide'))

def normalize_counter(x):
    '''normalize a counter object'''
    total = sum(x.values(), 0.0)
    for key in x:
        x[key] /= total
    return(x)