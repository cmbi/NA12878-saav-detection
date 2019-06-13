#!/usr/bin/env python3

import matplotlib, re, os,sys
from matplotlib_venn import venn3, venn3_unweighted
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sets import Set
from tqdm import tqdm
from collections import Counter

# Set the visualization settings (maybe need to be adjusted for saving figures to file)
matplotlib.rcParams['axes.titlesize'] = 'xx-large'
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

'''
this script will concatenate all the separate spectral match files generated from ionbot and process them to answer research questions

questions:
- horizontal/vertical coverage of the whole proteome?
- extent of multiple mapping, proportion of mapping that is exclusive to the long-read transcriptome
- are my peptides of interest present? How often does ionbot correctly predict an SNV?
'''
def concatenate_csvs(folder_path):
    directory= os.fsencode(csvpath)
    ionbotout=pd.DataFrame()
    for csvfile in os.listdir(directory):
        csvname=os.fsdecode(csvfile)
        if csvname.endswith('.csv'):
            temp=read_df_in_chunks(csvname, 1000)
            ionbotout=pd.concat([ionbotout,temp])
    return(ionbotout)

def read_df_in_chunks(filename, chunksize):
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv(filename, chunksize=chunksize) # chunksize represents number of rows read per chunk

    chunk_list = []  # append each chunk df here 

    # Each chunk is in df format
    for chunk in df_chunk:  
        # perform data filtering 
        #chunk_filter = chunk_preprocessing(chunk)

        # Once the data filtering is done, append the chunk to list
        chunk_list.append(chunk)

    # concat the list into dataframe 
    df_concat = pd.concat(chunk_list) # this is your final dataframe
    
    return(df_concat)

def chunk_preprocessing(df_chunk):
    new_chunk=df_chunk[(df_chunk['ri_126.1277']>0) & (df_chunk['q_value']<=0.01) & (df_chunk['DB']=='T')]
    new_chunk=new_chunk[['scan_id','charge','precursor_mass','matched_peptide','modifications','ionbot_psm_score','DB','unexpected_modification','ms2pip_pearsonr','proteins','num_unique_pep_ids']]
    return(new_chunk)

def import_coding_transcriptids(sources):
    '''
    input: paths of protein sequences from reference and cell-specific transcriptome translation
    output: ids of the combination
    
    '''
    transcript_ids=[]
    sources=['/data/data/genome_grch38/gencode.v29.pc_translations.fa','/data/data/spectra/dictionary/dumb_orfprediction_setA/flair.setA.final.pep']
    for f in sources:
        with open(f) as handle:
            for line in handle:
                if line.startswith('>'):
                    if ' ' in line.strip():
                        tid=line.split(' ')[0]
                        transcript_ids.append(tid[1:])
                    else:
                        transcript_ids.append(line.strip()[1:])
    return(transcript_ids)

def detected_proteins(ids,pco):
    proteins_covered=pco
    for idu in ids:
        proteins_covered[idu]+=1
    return(proteins_covered)

def import_cpdt(cpdt,wantFull):
    ''' read the cpdt files into a data structure
    this function can also handle the cpdt files generated with interesting_peptide_finder (only peptides with SNVs)
    {protein_ID:{pep1:0, pep2:0, pep3:0}}
    '''
    cpdt_pep={}
    full_seqs={}
    for cpf in source:
        with open(cpf) as c:
            for line in c:
                if line.startswith('>'):
                    key=line.strip()[1:]
                    key=get_id(key)
                    cpdt_pep[key]={}
                    full_seqs[key]=''
                elif 'PEPTIDE' in line:
                    lp=line.split('PEPTIDE ')[1]
                    lp=lp.split(':')[0]
                    cpdt_pep[key][lp]=0
                elif 'PEPTIDE' not in line:
                    full_seqs[key]=line.strip()
    if wantFull:
        return(cpdt_pep, full_seqs)
    return(cpdt_pep)

def import_gff(gfffile,isBed):
    '''use the gfffile to associate what proteins belong to which chromosome, in order to show the chromosomal distribution
    '''
    chromdict={}
    with open(gfffile) as handle:
        for line in handle:
            info=line.split('\t')
            if len(info)>5:
                if not isBed:
                    if 'transcript_type=protein_coding' in line:
                        tid=line.split('transcript_id=')[1]
                        tid=tid.split(';')[0]
                        chromdict[tid]=info[0]
                else:
                    chromdict[info[3]]=info[0]
    return(chromdict)

def find_chrom(prots,chromdict):
    for p in prots:
        if p in chromdict:
            return(chromdict[p])
    return("unknown")

def get_id(idstring):
    i=idstring
    if '|m.' in i:
        i=i.split('|m.')[0]
    elif 'ENSP' in i:
        i=i.split('|')[1]
    return(i)

def coverage_measure(cpdt_pep,full_seqs):
    #high_cov_vert={}
    high_cov_hor={}
    perc_cov_dist=[]
    vert_cov=[]
    for p,peps in cpdt_pep.items():
        seq=full_seqs[p]
        remains=seq
        count_pep=0
        for s,c in peps.items():
            if c>5: #what constitutes a "true" hit
                count_pep+=c
                if s in remains:
                    remains=remains.replace(s,'')
                else:
                    prefix=re.split('R|K',s)
                    for p in prefix:
                        if len(p)>3 and p in remains:
                            remains=remains.replace(p,'')
        perc_cov=float((len(seq)-len(remains))/len(seq))*100
        perc_cov_dist.append(perc_cov)
        vert_cov.append(count_pep)
        if perc_cov>50:
            high_cov_hor[p]=peps
        # if count_pep>100:
        #     high_cov_vert[p]=peps
    return(high_cov_hor,vert_cov,perc_cov_dist)

def bin_hits_by_source(scanid,ids,oro,ono,ob):
    '''sort peptide hits by their source dictionary'''
    ref_only=oro
    ont_only=ono
    both=ob
    ont=False
    ref=False
    for i in ids:
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
    else:
        raise Exception('Unexpected protein found')
    return(ref_only,ont_only,both)

def fill_cpdt(pep,mod,ids,old_cpdt_pep):
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
        if i in cpdt_pep.keys():
            ispeplist=cpdt_pep[i]#make copy to mutate as you iterate
            for p,ct in cpdt_pep[i].items():
                if p==pep:
                    ct+=1
                    ispeplist[p]=ct
                    found=True
            cpdt_pep[i]=ispeplist
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


def plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_combi):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution")
    sns.distplot(ibdf_refonly['ms2pip_pearsonr'], hist=False, label='Human reference only',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_ontonly['ms2pip_pearsonr'], hist=False, label='Transcriptome translation only',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_combi['ms2pip_pearsonr'], hist=False, label='Reference + transcriptome translation',axlabel='Pearson R Correlation')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr.png")
    return("Scores plot made")

def plot_source_piechart(ref_only,ont_only,both):
    '''this function will plot the source piechart of sources of the hits and save it to a pdf'''
    explode = (0.1, 0, 0)
    labels='Exclusively ONT transcriptome','Exclusively reference (gencode)', 'Both'
    plt.pie([len(ont_only),len(ref_only),len(both)],autopct='%1.1f%%', explode=explode,colors=['#de2d26','#3182bd','#756bb1'])
    plt.title('Peptide spectral hits by source',fontsize=35)
    plt.legend(labels,loc = (0.1,-0.25))
    plt.savefig("sources_spectral_hits.png") 
    return("saved to sources_spectral_hits")

def plot_chromosomal_dist(distr_list):
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    #horizontal coverage
    plt.hist(distr_list)
    plt.ylabel("# Peptides")
    plt.xlabel("Chomosomes")
    plt.savefig("chromosomal_distribution.png")
    return("plotted chromosomal distribution")

def plot_coverage_plots(cpdt_pep):
    '''this function will plot the graphs that correspond to the coverage of the proteome
    - vertical coverage
    - horizontal coverage
    - chromosome distribution
    '''
    high_cov_hor,cov_vert,perc_cov_dist=coverage_measure(cpdt_pep)
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    #horizontal coverage
    plt.hist(perc_cov_dist,bins=300)
    plt.xlim(1,100)
    plt.xlabel("% Coverage")
    plt.ylabel("# Proteins")
    plt.savefig("horizontal_coverage.png")
    plt.clf()
    #vertical coverage
    plt.hist(cov_vert,bins=500)
    plt.xlabel("# Proteins")
    plt.ylabel("# Peptides")
    plt.savefig("vertical_coverage.png")
    return("Plotted coverage")

def plot_mut(mutant_cpdtpep,cpdtpep):
    '''plot protein abundance vs number of detected mutant peptides'''
    prot_abundance=[]
    nr_mutant=[]
    mut_proteins_detected=set()
    num_peptides=0
    num_occurences=0
    for prot,peps in mutant_cpdtpep.items():
        if prot in cpdtpep:
            sum_mut=0 #total number of detected mutant peptides
            sum_nonmut=0
            for pep,ct in peps.items():
                sum_mut+=ct
                if ct>0:
                    mut_proteins_detected.add(prot)
                    num_occurences+=ct
                    num_peptides+=1
            if sum_mut>0: #only record the proteins with at least 1 detected mutation peptide
                nr_mutant.append(sum_mut)
                for pepc,ctc in cpdtpep[prot]:
                    sum_nonmut+=ctc
                prot_abundance.append(sum_mut+sum_nonmut)
    #make plot
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    ax=sns.scatterplot(prot_abundance,nr_mutant)
    ax.set(xlabel='Protein abundance',ylabel='Number mutant peptides detected')
    ax.savefig("mutant_abundance.png")
    return(len(mut_proteins_detected),num_peptides,num_occurences)

def make_report(hits_df):
    '''report general information collected about the hits, including quality control measures
    which includes:
    - how many matched peptides above threshold
    - how many unique peptides are found
    - how many unique proteins are found (and how many are in transcriptome v reference)
    - inventory PTMs (including how many SNVs there are)
    - how many hits were removed with the quality filter
    - how many hits passing the quality threshold were not included because not matching with a theoretical peptide
    '''
    uniquepeptides=hits_df['peptide'].nunique()
    return("report written")

    
def main(directory_ontonly, directory_refonly, directory_combination, cpdtfile,cpdtfile_mut,gfffile,bedfile):
    '''this is the main function that will iterate over the giant pandas df and perform all analyses and make all figures
    this is written in a way that iteration should only be done once.
    '''
    #imports
    print("importing data")
    ibdf_ontonly=concatenate_csvs(directory_ontonly)
    ibdf_refonly=concatenate_csvs(directory_refonly)
    ibdf_combi=concatenate_csvs(directory_combination)

    #qc function
    plot_scores(ibdf_ontonly,ibdf_refonly,ibdf_combi)

    #filter badly scoring hits
    ibdf_ontonly = chunk_preprocessing(ibdf_ontonly)
    ibdf_refonly = chunk_preprocessing(ibdf_refonly)
    ibdf_combi = chunk_preprocessing(ibdf_combi)

    #import insilico digest info
    cpdt_pep,full_seqs=import_cpdt(cpdtfile,True) #import cpdt will all peptides (cat gencode and flair beforehand)
    mut_cpdt_pep=import_cpdt(cpdtfile_mut,False) #import the cpdt file with all snv peptides
    chromdict_ref=import_gff(gfffile,False) #import gff3 file annotations from gencode
    chromdict_ont=import_gff(bedfile,True) #import bed file annotations from ont (converted from psl)
    chromdict={**chromdict_ont,**chromdict_ref} #combine the 2 dictionaries

    
    #initialize data structures to collect information
    print("initializing analyses...")
    proteins_covered=Counter() #proteins detected
    mutated=set() #all protiens that matched to a predicted mutated peptide
    ref_only=set() #scan ids in the reference set
    ont_only=set() #scan ids in the ont set
    both=set() #scan ids that matched to both ref and ont proteins
    chrom_dist=[] # 1 chromosome location per scan id
    hits_missed=0
    hit_mut=0
    

    #iterate to fill the data structures
    print("Analyzing data...")
    for row in tqdm(ibdf_combi.iterrows()):
        scanid=row[1][0]
        mod=str(row[1][7])
        pep=row[1][3]
        if '||' in row[1][9]:
            ids=row[1][9].split('||')
        else:
            ids=[row[1][9]]
        proteins_covered=detected_proteins(ids,proteins_covered) #what proteins from the proteome are covered and in what amounts
        ref_only,ont_only,both=bin_hits_by_source(scanid,ids,ref_only,ont_only,both) #what dictionaries do the hits come from
        cpdt_pep,notfound=fill_cpdt(pep,mod,ids,cpdt_pep) #what mutant peptides are detected, how many, and what proteins they come from
        hits_missed+=notfound
        chrom_dist.append(find_chrom(ids,chromdict)) #which chromosome does the 

        if '->' in mod: #if ib detects a mutated peptide
            hit_mut+=1
            mut_cpdt_pep=fill_cpdt(pep,mod,ids,mut_cpdt_pep)
            for i in ids:
                mutated.add(get_id(i))

    #create the figures
    print("number of hits with detected mutation = " +str(len(hit_mut))+ " matched to "+str((mutated))+ " proteins.")
    print("number of hits that were not counted because they were not predicted by in silico digest: "+str(hits_missed))
    mut_proteins_detected,mut_peptides,mut_occurences=plot_mut(mut_cpdt_pep,cpdt_pep)
    print("Total of "+str(mut_occurences)+" occurances of "+str(mut_peptides)+" peptides from "+str(mut_proteins_detected)+" proteins were detected")
    plot_coverage_plots(cpdt_pep)
    plot_source_piechart(ref_only,ont_only,both)
    plot_chromosomal_dist(chrom_dist)
    return("Finished")
    
    
main(*sys.argv[1:])
    
    
    
    
    