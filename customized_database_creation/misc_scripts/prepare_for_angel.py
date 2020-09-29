#!/usr/bin/env python2

import re

def check_actg(ceoutput):
    """takes the combine_exons output and checks that the transcripts that are there only have ACTG
    """
    wrong=0
    all=0
    with open(ceoutput) as handle:
        for line in handle:
            if not line.startswith(">"):
                all+=1
                if len(re.findall('[^A|^T|^G|^C]',line.strip()))>0:
                    wrong+=1
    return str(wrong) + ' out of '+str(all)+ ' have a non ATCG character'

def remove_non_actg(ceoutput,outfile):
    """this function removes the 2 sequences that have non-actg characters in them
    """
    f=open(outfile,'w')
    seq1=''
    with open(ceoutput) as handle:
        for line in handle:
            if line.startswith(">"):
                header=line
            elif len(re.findall('[^A|^T|^G|^C]',line.strip()))==0:
                f.writelines(header+line)
    f.close()
    return "reached end"

print remove_non_actg(sys.argv[1],sys.argv[2])