#!/usr/bin/env python3

'''
this script will use as input:
    - the report output from phase_exons
    - CDS from both ref and ont
    - file with identifiers of the proteins and peptides from fdr_reestimation

will return the SNP information (first file) from those variants that were found in the peptides (last file)

1. parse CDS (only include observed) from ANGEL output (suffix "final.cds")
2. add column to the report with protein ids
3. get sequence
'''

def read_filter_cds(cdsfile,prots):
    cds=pd.DataFrame(columns=["protein","cds"])
    with open(cdsfile) as handle:
        for line in handle:
            if line.startswith('>'):
                add=False
                header=line.split('|m.')[0]
                if header[1:] in prots: #check if this CDS is important to keep
                    header=header[1:]
                    add=True
            else:
                seq=line.strip()
                if add:
                    cds=cds.append({"protein":header,"cds":seq})
    return(cds)



def map_seq(reportdf):
    '''
    separate the list of proteins, so that it is 1 protein, 1 peptide format
    '''

def gather_prots(list_prots):
    out=set()
    for p in list_prots:
        if '||' in p:
            s=set(p.split('||'))
            out=out.union(s)
        else:
            out.add(p)
    return out

def protein_filter(prots):
    o=set()
    for p in prots:
        o.add(get_id(p))
    return(o)

def get_id(idstring):
    i=idstring
    if '|m.' in i:
        i=i.split('|m.')[0]
    elif 'ENSP' in i:
        i=i.split('|')[1]
    return(i)

def main():
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--var', help='', required=True)
    parser.add_argument('--cdsont', help='', required=True)
    parser.add_argument('--cdsref', help='', required=True)
    parser.add_argument('--observed', help='', required=True)
    args=parser.parse_args()

    #read in the observed
    obs=pd.read_csv(args.observed,sep='\t') #proteins \t peptide
    # obs=obs.str.split(',',expand=True).stack()
    prots_unfiltered=obs["proteins"].unique()
    prots=gather_prots(prots_unfiltered)
    cds_ont=read_filter_cds(args.cdsont,prots) #protein \t cds
    cds_ref=read_filter_cds(args.cdsref,prots)
    report=pd.read_csv(args.var,sep='\t',columns=["transcript","local_pos","chrom_pos"]) #transcript \t var positions on transcript \t chr|chrom_positions
    transcripts=protein_filter(prots)
    report=report.loc[report["transcript"].isin(transcripts)]
    
    