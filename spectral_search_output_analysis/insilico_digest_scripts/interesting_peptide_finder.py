#!/usr/bin/env python2

import re, collections, sys
from sets import Set
import multiprocessing as mp
import argparse

def insilico_digest_diff(cpdtref,cpdtcustom,cpu):
    '''this script will find peptides that vary between two cpdt files
    prints cpdt files: one for all different peptides and one with only snv peptides
    '''
    print('Reading in files...')
    ref=read_csv(cpdtref)
    custom=read_csv(cpdtcustom)
    print('Searching for variant peptides...')
    pool=mp.Pool(processes=int(cpu),initializer=child_initialize, initargs=(ref,custom))
    results=[pool.apply_async(snvfinder,args=(pid,peplist,)) for pid,peplist in ref.iteritems()]
    output = [p.get() for p in results] #list of tuplies/list of dictionaries
    print('Gathering output...')
    return combine_output(output)

def child_initialize(_ref,_custom):
     global ref, custom
     ref = _ref
     custom = _custom
    
def snvfinder(pid,peplist):
    idho=pid+'_h1'
    idhz=pid+'_h0'
    keys=[]
    newsnv={}
    newall={}
    snvcounterpart={}
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
                intersect=c_peplist.intersection(r_peplist)
                if len(intersect)!=0: # should have some peptides in common
                    for pep in dif:
                        if isQualified(pep,ref): #need to check if variant peptide is not already present in the reference dictionary
                            prob=pepdict[pep] #fetch start position (not probability)
                            # if float(prob)>=0.05: #higher than cutoff probability
                            isSNV,counterpart=determine_snv(pep,r_peplist)
                            if isSNV: #counterpart not filtered (optionally can filter that it does not show up in another gene)
                                if key not in newsnv:
                                    newsnv[key]={}
                                    snvcounterpart[key]={}
                                newsnv[key][pep]=prob
                                snvcounterpart[key][counterpart]=peplist[counterpart]
                            if key in newall:
                                newall[key][pep]=prob
                            else:
                                newall[key]={pep:prob}
    return [newsnv,newall,snvcounterpart]

def combine_output(process_output):
    '''
    takes the multiprocessing output and concatenates the appropriate libraries
    '''
    newsnv,newall,snvcounterpart={},{},{}
    for l in process_output:
        newsnv=merge_two_dicts(newsnv, l[0])
        newall=merge_two_dicts(newall,l[1])
        snvcounterpart=merge_two_dicts(snvcounterpart,l[2])
    return newsnv,newall,snvcounterpart

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def get_id(idstring):
    if '|m.' in idstring:
        outstring=idstring.split('|m.')[0]
    elif 'ENSP' in idstring:
        tid=idstring.split('|')[1]
        if 'Random' in idstring:
            prefix=idstring.split('_')[0]
            outstring=prefix+'_'+tid
        else:
            outstring=tid
    else:
        outstring=idstring.strip()
    return(outstring)

def isQualified(peptide,reference_pepdict,reference=False):
    count=0
    for pid,peplist in reference_pepdict.iteritems():
        if peptide in peplist:
            count+=1
            if reference and count>1:
                return False
            elif not reference:
                return False
    return True

def determine_snv(peptide,plist):
    ''' checks whether the peptide in question differs from a member in the list by exactly 1 amino acid
    input: a peptide and a list of peptides
    output: boolean
    '''
    mismatch=0
    il=False #TO IMPLEMENT: get statistics about how many times we have an IL substitution as the only substitution
    for pep in plist:
        if len(pep)==len(peptide) and pep!=peptide:
            for idx,aa in enumerate(pep):
                if aa!=peptide[idx]:
                    situation1= aa=='I' and peptide[idx]=='L'
                    situation2= aa=='L' and peptide[idx]=='I'
                    if not situation1 and not situation2:
                        mismatch+=1
                    else:
                        il=True
            if mismatch==1:
                return (True,pep)
    return (False,'')

def read_cpdt(cpdt):
    '''{id:{peptide:probability}'''
    cpdt_pep={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=get_id(line)
                cpdt_pep[key]={}
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(': ')
                cpdt_pep[key][lp[0]]=lp[1].strip()
    return cpdt_pep

def read_csv(csvfile):
    cpdt_pep={}
    with open(csvfile) as c:
        for line in c:
            if 'protein' not in line:
                info=line.split(',')
                key=get_id(info[0])
                if key not in cpdt_pep:
                    cpdt_pep[key]={}
                cpdt_pep[key][info[1]]=info[2].strip()
    return cpdt_pep

def write_csv(d,outfile):
    '''
    input format: {id:{peptide:probability}}
    output: new csv
    '''
    f=open(outfile,'w')
    f.writelines("{},{},{}".format('protein', 'matched_peptide', 'start')+'\n')
    for id,peplist in d.iteritems():
        # f.writelines(id+'\n')
        for pep,prob in peplist.iteritems():
            f.writelines(id+','+pep+','+prob+'\n')
    return 'new CPDT file written to '+outfile

def main():
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--cpu',help='CPU number',required=True)
    parser.add_argument('--ref', help='Reference csv file', required=True)
    parser.add_argument('--var', help='Variant containing csv file', required=True)
    parser.add_argument('--outall', help='Output file all differing', required=False)
    parser.add_argument('--outsnp', help='Output file snv differing', required=False)
    parser.add_argument('--outctp', help='Output file reference counterpart (to snv)', required=False)
    args=vars(parser.parse_args())
    snvdiffpeps,alldiffpeps,snvcounterpart=insilico_digest_diff(args['ref'],args['var'],args['cpu'])
    print('Writing to files...')
    if args['outall']:
        write_csv(alldiffpeps,args['outall'])
    if args['outsnp']:
        write_csv(snvdiffpeps,args['outsnp'])
    if args['outctp']:
        write_csv(snvcounterpart,args['outctp'])

if __name__ == "__main__":
    main()