#!/usr/bin/env python3

import random, argparse, sys
from pyteomics.parser import cleave

"""
prepare an appropriate decoy set for re-estimation of FDR for variant peptides
1. Find out the positions of the variant/counterpart peptides, as well as position on the reversed sequences
2. Reverse the sequences in the databases
3. Infer the CPDT peptide from the sequence+position on the reversed sequence
"""

def import_cpdt(cpdt):
    ''' read the cpdt files into a data structure
    this function can also handle the cpdt files generated with interesting_peptide_finder (only peptides with SNVs)
    {protein_ID:[pep1, pep2, pep3]} 
    '''
    cpdt_pep={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=get_id(key)
                cpdt_pep[key]=[]
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(':')
                cpdt_pep[key].append(lp[0])
    return(cpdt_pep)

def import_seq(cpdt):
    full_seqs={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=get_id(key)
                full_seqs[key]=''
            elif 'PEPTIDE' not in line:
                full_seqs[key]=line.strip()
    return(full_seqs)

def get_id(idstring):
    i=idstring
    if '|m.' in i:
        i=i.split('|m.')[0]
    elif 'ENSP' in i:
        i=i.split('|')[1]
    return(i)

def write_cpdt(d,outfile):
    '''
    input format: {id:{peptide:probability}}
    output: new cpdt file
    '''
    f=open(outfile,'w')
    for id,peplist in d.items():
        f.writelines('>'+id+'\n')
        for pep in peplist:
            f.writelines('\t'+'PEPTIDE '+pep+': \n')
    return('new CPDT file written to '+outfile)

def get_rev_peps(cpdt_var,cpdt_ctp,full_seqs,missed_cleavages):
    '''
    find the variant on the peptide and get the position on the peptide
        given:
        - {protein_id:[var_peps]}
        find:
        - {var_pep:pep_position} with counterpart info
        - {var_pep:start_position_var_pep} with sequence info
    then find the position of the peptide on the sequence + position on peptide = position of variant on sequence
    then len(sequence)-position of variant on sequence= position of variant on reversed sequence
    then find the peptide(s) in the reverse digest that it would be in
    '''
    new_cpdt={}
    for prot, plist in cpdt_var.items():
        plist_ctp=cpdt_ctp[prot]
        seq=full_seqs[prot]
        rev_seq=seq[::-1]
        prelim_peps=[]
        for pep in plist:
            pos_pep=var_pos_on_pep(pep,plist_ctp)
            pos_seq=pep_pos_on_seq(pep,seq)
            pos_var=pos_pep+pos_seq
            pos_var_rev=len(seq)-pos_var
            seq_begin,seq_end=rev_seq[:pos_var_rev],rev_seq[pos_var_rev:]
            if seq_begin.rfind('K')>seq_begin.rfind('R'):
                start_revpep=seq_begin.rfind('K')
            elif seq_begin.rfind('R')>seq_begin.rfind('K'):
                start_revpep=seq_begin.rfind('R')
            else:
                start_revpep=0
            first_half=seq_begin[start_revpep:]
            if seq_end.find('K')>seq_end.find('R'):
                end_revpep=seq_end.find('R')
            elif seq_end.find('R')>seq_end.find('K'):
                end_revpep=seq_end.find('K')
            else:
                end_revpep=len(seq)
            second_half=seq_end[:end_revpep]
            target_pep=first_half+second_half
            prelim_peps.append(target_pep[1:])
        #find the actual pep from a digest that corresponds to the target pep
        cleaved=cleave(rev_seq,'[KR]',int(missed_cleavages),4)
        final_peps=match_prelim_peps(prelim_peps,cleaved)   
        new_cpdt[prot]=final_peps
    return new_cpdt

def fetch_rev_counterparts(new_cpdt,ref_seq,missed_cleavages):
    '''
    for the variant-free library, no variant position is present.
    however, the reference counterparts of the just-found reverse variant peptides can be used as the decoy for this set

    input: the cpdt file output from get_rev_peps (containing all the peptides used as decoy for the variant containing set) and the reference sequences in a dictionary
    '''
    var_free_decoys={}
    cpdt_ref=rev_theoretical_digest(ref_seq,missed_cleavages) #theoretical digest of the reversed reference proteins
    for prot,peplist in new_cpdt.items():
        new_peplist=[]
        if '_h' in prot:
            prot=prot.split('_h')[0]
        for pep in peplist:
            print(pep,cpdt_ref[prot])
            sys.exit()
            isSNV,counterpart=determine_snv(pep,cpdt_ref[prot])
            if isSNV:
                new_peplist.append(counterpart)
            elif pep in cpdt_ref[prot]:
                new_peplist.append(pep)
    return var_free_decoys

def rev_theoretical_digest(seq_dict,missed_cleavages):
    digest={}
    for prot,seq in seq_dict.items():
        seq_cut = cleave(seq, '[KR]', int(missed_cleavages),4)
        digest[prot]=seq_cut
    return digest

def determine_snv(peptide,plist):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: boolean
    '''
    mismatch=0
    for pep in plist:
        if len(pep)==len(peptide) and pep!=peptide:
            for idx,aa in enumerate(pep):
                if aa!=peptide[idx]:
                    situation1= aa=='I' and peptide[idx]=='L'
                    situation2= aa=='L' and peptide[idx]=='I'
                    if not situation1 and not situation2:
                        mismatch+=1
            if mismatch==1:
                return(True,pep)
    return(False,'')

def var_pos_on_pep(pep,counterpart):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: position
    '''
    for c in counterpart:
        if len(pep)==len(c): #this should always be true
            original=''
            sub=''
            mismatch=0
            position=0
            for idx,aa in enumerate(pep):
                if aa!=c[idx]:
                    original=aa
                    sub=c[idx]
                    situation1= original=='I' and sub=='L'
                    situation2= original=='L' and sub=='I'
                    if not situation1 and not situation2:
                        mismatch+=1
                        position=idx
            if mismatch==1:
                return position
    return(0)

def pep_pos_on_seq(peptide,seq):
    pieces=seq.split(peptide)
    return(len(pieces[0]))

def match_prelim_peps(list_prelim, list_digest):
    list_final=[]
    for pep in set(list_prelim):
        for pepe in list_digest:
            if pep in pepe:
                list_final.append(pepe)
    return(list_final)

def main():
    parser= argparse.ArgumentParser(description="Variant peptide decoy generator")
    parser.add_argument('--vcpdt',type=str, help='SAAV CP-DT file with fasta header included')
    parser.add_argument('--ccpdt',type=str, help='reference counterpart CP-DT file with fasta header included')
    # parser.add_argument('--refseqcpdt',type=str, help="a cp-dt file with header and the sequence included (no variants)")
    parser.add_argument('--varseqcpdt',type=str, help="a cp-dt file with header and the sequence included (with variants)")
    # parser.add_argument('--mc', help='Missed cleavages')
    args = vars(parser.parse_args()) 
    #read in the cpdts/seqs
    cpdt_var=import_cpdt(args["vcpdt"])
    cpdt_ctp=import_cpdt(args["ccpdt"])
    # ref_seqs=import_seq(args["refseqcpdt"])
    var_seqs=import_seq(args["varseqcpdt"])
    #create the decoy lists
    decoys_var_cont=get_rev_peps(cpdt_var,cpdt_ctp,var_seqs,1)
    # decoys_var_free=fetch_rev_counterparts(decoys_var_cont,ref_seqs,1) #this doesn't work
    #write to files
    write_cpdt(decoys_var_cont,'var_cont_decoys.cpdt')
    # write_cpdt(decoys_var_free,'var_free_decoys.cpdt')

if __name__ == '__main__':
	main()
