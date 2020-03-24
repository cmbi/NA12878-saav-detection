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
    return(pd.DataFrame({'sub':all_list}))

def bl_62():
    return({
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4})

# def counter_to_df(allc, observed):
#     '''create the df that will be used in the heatmap figure, use normalization (min-max)'''
#     df=get_normalized_matrix(allc,observed)
#     return(df,dfall)

def determine_snv(peptide,counterpart):
    ''' extract the single amino acid substitution from the 
    input: a peptide and a list of peptides
    output: tuple
    ***OUTDATED, NO LONGER IN USE***
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

def categorize_mods(mod,pred_aa_sub):
    saav=False
    if pred_aa_sub!='':
        saav=True
    if mod!='':
        if saav:
            return('SAAV + mod')
        else:
            return(mod)
    return('SAAV' if saav else 'none')

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

def normalize_counter(x):
    '''normalize a counter object'''
    total = sum(x.values(), 0.0)
    for key in x:
        x[key] = (x[key]/total) * 100
    return(x)

def match_var_nonvar(df_var,df_nonvar,var_pep_df):
    pep_ct_var=df_var['variant_peptide'].value_counts().reset_index()
    pep_ct_var.columns=['variant_peptide','count_var']
    var_pep_df=var_pep_df.merge(pep_ct_var,on='variant_peptide',how='left')
    pep_ct_refctp=df_nonvar['ref_counterpart'].value_counts().reset_index()
    pep_ct_refctp.columns=['ref_counterpart','count_refctp']
    var_pep_df=var_pep_df.merge(pep_ct_refctp,on='ref_counterpart',how='left').fillna(0)
    return(var_pep_df[var_pep_df['count_var']!=0].drop(columns=['id','haplotype']).drop_duplicates())

def sub_conversion(list_aa):
    if len(list_aa)==2:
        codes={'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K','Trp': 'W', 'Asn': 'N', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Ala': 'A','Gly': 'G', 'Ile': 'X', 'Leu': 'X', 'Xle': 'X', 'His': 'H', 'Arg': 'R', 'Met': 'M','Val': 'V', 'Glu': 'E', 'Tyr': 'Y'}
        new_list=[]
        for l in list_aa:
            new_list.append(codes[l])
        return(','.join(new_list[::-1])) #can also be a tuple later if desired
    return('')