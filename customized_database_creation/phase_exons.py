#!/usr/bin/env python2

import re, collections, sys, vcf
from sets import Set
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse

def phase_exons(fasta_exons,vcf_gz,outputfile):
    '''
    this function goes through the fasta files of the exons (made from getfasta of the gff3 + ref) and replaces variants
    if the exon contains at least 1 heterozygous variant, a 0 and 1 version of the exon is made
    includes correction for coordinates- vcf uses 1-based coordinates, bed uses 0
    the fasta has a the sequences on one line, no breaks
    
    input: fasta file of exons created from getfasta, vcf file of heterozygous variants excluding header
    output: new fasta files of exons, with 0 and 1 haplotypes for the exons that include one or more heterozygous variants
        header will include information about homozygous and heterozygous variants

    also print out report of what exons have what variants to make it easier to trace. This report will be tab delimited and consist of:
        - chr start-stop in the previous format
        - chromosomal position of variant
        - ref base
    '''
    f=open(outputfile,'w')
    print "loading variants..."
    variants=vcf.Reader(filename=vcf_gz)
    #variants=pd.read_csv(vcf, sep="\t", names=["chrom","pos","id","ref","alt","qual","filter","info","format","na12878"])
    print "replacing exon sequences..."
    f.writelines('\t'.join(['chromosome','start','end','haplotype','pos_het','het_org','pos_hom','hom_org','sequence'])+'\n')
    with open(fasta_exons) as fe:
        header=''
        entry_0=''
        entry_1=''
        bordervar=0
        varall=0
        het=False
        var=False
        keys=Set()
        for line in fe:
            if line.startswith('>'):
                header=line
            else:
                if header[1:].strip() not in keys: #prevent duplicates
                    keys.add(header[1:].strip())
                    #parse header
                    head=re.split(':|-',header[1:].strip())
                    chr,start,end=head[0],int(head[1]),int(head[2])
                    if chr!='chrY' and chr!='chrM': #the sample is a female, and there are no mitochondial variants in the vcf file
                        #retrieve appropriate vcf lines
                        vcf_fetch=variants.fetch(chr,start-1,end) #since the vcf file is 1-based
                        var_pos_end=0
                        seq=line.strip()
                        entry_0,entry_1='',''
                        het=False
                        var=False
                        hom_variants=[] #homozygous variants
                        hom_variants_org=[] #homozygous variant origins
                        het_variants=[] #heterozygous variants
                        het_variants_org=[] #heterozygous variant origins
                        for vari in vcf_fetch:
                            var_pos=max(0,vari.POS-start-1) #start position on the sequence string, don't allow neagtive
                            entry_0+=seq[var_pos_end:var_pos] #add portion between old end position and new start position
                            entry_1+=seq[var_pos_end:var_pos]
                            var_pos_end=var_pos+len(vari.REF) #assign new ending position
                            if var_pos_end>len(seq):
                                var_pos_end=len(seq)
                            #add a check that the base(s) that I am replacing are what they should be
                            ref_allele=seq[var_pos:var_pos_end]
                            if ref_allele==vari.REF: #only if the sequence that is being replaced matches what is written in the vcf file. this discludes all border
                                varall+=1
                                if vari.genotype('NA12878')['GT']=="0|1":
                                    het=True
                                    het_variants.append(str(var_pos))
                                    het_variants_org.append(chr+'|'+str(vari.POS)) 
                                    entry_0+=vari.REF
                                    entry_1+=str(vari.ALT[0])
                                elif vari.genotype('NA12878')['GT']=="1|0":
                                    het=True
                                    het_variants.append(str(var_pos))
                                    het_variants_org.append(chr+'|'+str(vari.POS))
                                    entry_0+=str(vari.ALT[0])
                                    entry_1+=vari.REF
                                elif vari.genotype('NA12878')['GT']=="1|2":
                                    het=True
                                    het_variants.append(str(var_pos))
                                    het_variants_org.append(chr+'|'+str(vari.POS))
                                    entry_0+=str(vari.ALT[0])
                                    entry_1+=str(vari.ALT[1])
                                elif vari.genotype('NA12878')['GT']=="2|1":
                                    het=True
                                    het_variants.append(str(var_pos))
                                    het_variants_org.append(chr+'|'+str(vari.POS))
                                    entry_0+=str(vari.ALT[1])
                                    entry_1+=str(vari.ALT[0])
                                elif vari.genotype('NA12878')['GT']=="1|1":
                                    var=True
                                    hom_variants.append(str(var_pos))
                                    hom_variants_org.append(chr+'|'+str(vari.POS))
                                    entry_0+=str(vari.ALT[0])
                                    entry_1+=str(vari.ALT[0])
                                else:
                                    print 'unexpected input from the final column of the vcf file '+ vari.genotype('NA12878')['GT']
                                    break
                            else:
                                #add back what was already there- no replacement (in case of border variants
                                entry_0+=ref_allele
                                entry_1+=ref_allele
                                #count the border variants
                                bordervar+=1
                        entry_0+=seq[var_pos_end:end] #add last chunk of sequence
                        entry_1+=seq[var_pos_end:end]
                        if het: #if contains heterozygous, make 2 seperate sequences
                            line_zero=[head[0],head[1],head[2],str(0),','.join(het_variants),','.join(het_variants_org)]
                            line_one=[head[0],head[1],head[2],str(1),','.join(het_variants),','.join(het_variants_org)]
                            if var: #if also contains homozygous, add those variants too
                                line_zero.extend([','.join(hom_variants),','.join(hom_variants_org)])
                                line_one.extend([','.join(hom_variants),','.join(hom_variants_org)])
                            else:
                                line_zero.extend([str(None),str(None)])
                                line_one.extend([str(None),str(None)])
                            line_zero.append(entry_0)
                            line_one.append(entry_1)
                            # f.writelines(header_0+'\n'+entry_0+'\n'+header_1+'\n'+entry_1+'\n')  
                            f.writelines('\t'.join(line_zero)+'\n'+'\t'.join(line_one)+'\n')
                        elif var: #if only contains homozygous variants
                            # f.writelines(header.strip()+' pos_hom:'+','.join(v)+'\n'+entry_0+'\n')
                            f.writelines('\t'.join([head[0],head[1],head[2],str(None),str(None),str(None),','.join(hom_variants),','.join(hom_variants_org),entry_0])+'\n')
                        else: #no variants in the exon
                            # f.writelines(header+line)
                            f.writelines('\t'.join([head[0],head[1],head[2],str(None),str(None),str(None),str(None),str(None),line.strip()])+'\n')
    f.close()
    return 'number total variants replaced = '+str(varall)+'. total border variants that were skipped = '+str(bordervar)

# def get_variant_sequences(chromosome,start,end,seq):
#     vcf_fetch=variants.fetch(chromosome,start-1,end) #since the vcf file is 1-based
#     var_pos_end=0
#     entry_0,entry_1='',''
#     hom_variants=[] #homozygous variants
#     hom_variants_org=[] #homozygous variant origins
#     het_variants=[] #heterozygous variants
#     het_variants_org=[] #heterozygous variant origins
#     for vari in vcf_fetch:
#         var_pos=max(0,vari.POS-start-1) #start position on the sequence string, don't allow neagtive
#         entry_0+=seq[var_pos_end:var_pos] #add portion between old end position and new start position
#         entry_1+=seq[var_pos_end:var_pos]
#         var_pos_end=var_pos+len(vari.REF) #assign new ending position
#         if var_pos_end>len(seq):
#             var_pos_end=len(seq)
#         #add a check that the base(s) that I am replacing are what they should be
#         ref_allele=seq[var_pos:var_pos_end]
#         if ref_allele==vari.REF: #only if the sequence that is being replaced matches what is written in the vcf file. this discludes all border
#             varall+=1
#             if vari.genotype('NA12878')['GT']=="0|1":
#                 het_variants.append(str(var_pos))
#                 het_variants_org.append(chr+'|'+str(vari.POS)) 
#                 entry_0+=vari.REF
#                 entry_1+=str(vari.ALT[0])
#             elif vari.genotype('NA12878')['GT']=="1|0":
#                 het_variants.append(str(var_pos))
#                 het_variants_org.append(chr+'|'+str(vari.POS))
#                 entry_0+=str(vari.ALT[0])
#                 entry_1+=vari.REF
#             elif vari.genotype('NA12878')['GT']=="1|2":
#                 het_variants.append(str(var_pos))
#                 het_variants_org.append(chr+'|'+str(vari.POS))
#                 entry_0+=str(vari.ALT[0])
#                 entry_1+=str(vari.ALT[1])
#             elif vari.genotype('NA12878')['GT']=="2|1":
#                 het_variants.append(str(var_pos))
#                 het_variants_org.append(chr+'|'+str(vari.POS))
#                 entry_0+=str(vari.ALT[1])
#                 entry_1+=str(vari.ALT[0])
#             elif vari.genotype('NA12878')['GT']=="1|1":
#                 hom_variants.append(str(var_pos))
#                 hom_variants_org.append(chr+'|'+str(vari.POS))
#                 entry_0+=str(vari.ALT[0])
#                 entry_1+=str(vari.ALT[0])
#             else:
#                 print 'unexpected input from the final column of the vcf file '+ vari.genotype('NA12878')['GT']
#                 break
#         else:
#             #add back what was already there- no replacement (in case of border variants
#             entry_0+=ref_allele
#             entry_1+=ref_allele
#             #count the border variants
#             bordervar+=1
#     entry_0+=seq[var_pos_end:end] #add last chunk of sequence
#     entry_1+=seq[var_pos_end:end]
#     return([','.join(hom_variants),','.join(hom_variants_org),','.join(het_variants),','.join(het_variants_org,entry_0,entry_1])

def main():
    parser = argparse.ArgumentParser(description='phase exons')
    parser.add_argument('--cds', help='CDS fasta sequences acquired by getfasta in bedtools', required=True)
    parser.add_argument('--vcf', help='VCF in gzipped format (vcf.gz)', required=True)
    parser.add_argument('--out', help='Output file for phased exons', required=True)
    args=vars(parser.parse_args())
    # cds=pd.read_csv(args['cds'],delim_whitespace=True,names=['chromosome','sequence']).drop_duplicates()
    # info=cds['chromosome'].str.split(':',expand=True)
    # cds['chromosome']=info[0]
    # positions=info[1].str.split('-')
    # cds['start']=positions[0]
    # cds['end']=positions[1]
    # cds['pos_het']=cds.apply(lambda x: get_variant_sequences(x['chromosome'],x['start'],x['end'],x['sequence']) ,axis=1)
    # df[['pos_het','het_org','pos_hom','hom_org','sequence_0','sequence_1']]=pd.DataFrame(df['pos_het'].values.tolist(),index=df.index)
    phase_exons(args['cds'],args['vcf'],args['out'])

if __name__ == "__main__":
    main()
