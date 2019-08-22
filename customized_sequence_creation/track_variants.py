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

def read_filter_cds(cdsfile,observed):
    prots_unfiltered=observed["proteins"].unique()
    prots=gather_prots(prots_unfiltered)
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

def main():
    parser = argparse.ArgumentParser(description='track variants')
    parser.add_argument('--exonvar', help='', required=True)
    parser.add_argument('--cdsont', help='', required=True)
    parser.add_argument('--cdsref', help='', required=True)
    parser.add_argument('--observed', help='', required=True)
    args=parser.parse_args()

    #read in the observed
    obs=pd.read_csv(args.observed,sep='\t') #proteins \t peptide
    obs=obs.str.split(',',expand=True).stack()
    cds_ont=read_filter_cds(args.cdsont,obs) #protein \t cds
    cds_ref=read_filter_cds(args.cdsref,obs)
    #need protein to exon mapping- with ref there is 1 line per exon mapping. with ont there is bed12
    report=read_csv #chr:pos-pos \t chromosomal variant position \t reference allele \t position on the exon sequence \t exon sequence
