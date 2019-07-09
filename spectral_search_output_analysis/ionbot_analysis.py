#!/usr/bin/env python3

import matplotlib, re, os,sys
from matplotlib_venn import venn2,venn2_unweighted,venn3, venn3_unweighted
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from collections import Counter
import json

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
def concatenate_csvs(csvpath):
    directory= os.fsencode(csvpath)
    ionbotout=pd.DataFrame()
    for csvfile in os.listdir(directory):
        csvname=os.fsdecode(csvfile)
        if csvname.endswith('.csv'):
            temp=read_df_in_chunks(os.path.join(csvpath,csvname), 1000)
            ionbotout=pd.concat([ionbotout,temp])
    return(ionbotout)

def read_df_in_chunks(directory, chunksize):
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv(directory, chunksize=chunksize) # chunksize represents number of rows read per chunk
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

def import_cpdt(cpdt):
    ''' read the cpdt files into a data structure
    this function can also handle the cpdt files generated with interesting_peptide_finder (only peptides with SNVs)
    {protein_ID:{pep1:0, pep2:0, pep3:0}}
    '''
    cpdt_pep={}
    full_seqs={}
    with open(cpdt) as c:
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
    return(cpdt_pep, full_seqs)

def import_cpdt_simple(cpdt):
    cpdt_pep={}
    with open(cpdt) as c:
        for line in c:
            if line.startswith('>'):
                key=line.strip()[1:]
                key=get_id(key)
                cpdt_pep[key]=set()
            elif 'PEPTIDE' in line:
                lp=line.split('PEPTIDE ')[1]
                lp=lp.split(':')[0]
                cpdt_pep[key].add(lp)
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
        p=get_id(p)
        if '_h' in p:
            p=p.split('_h')[0]
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
        vert_cov.append(float(count_pep/len(peps.keys())))
        if perc_cov>50:
            high_cov_hor[p]=peps
        # if count_pep>100:
        #     high_cov_vert[p]=peps
    return(high_cov_hor,vert_cov,perc_cov_dist)

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

def calculate_support():
    '''look into the support for proteins
    how many proteins have 1,2,3... peptides supporting their existence? unambiguously?'''
    return(0)

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
                break
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
    sns.distplot(ibdf_refonly['ionbot_psm_score'], hist=False, label='Human reference only',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_ontonly['ionbot_psm_score'], hist=False, label='Transcriptome translation only',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_combi['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation',axlabel='Pearson R Correlation')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr_3source.png")
    plt.clf()
    return("Scores plot made")

def plot_scores_pg(ibdf_combi,ibdf_combi_pg):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_combi['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation (no variants)',axlabel='Pearson R Correlation')
    sns.distplot(ibdf_combi_pg['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation (with variants)',axlabel='Pearson R Correlation')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr_pgvsopenmut.png")
    plt.clf()
    return("Scores plot made")

def plot_source_piechart(ref_only,ont_only,both,figname,isOpenmut):
    '''this function will plot the source piechart of sources of the hits and save it to a pdf'''
    plt.figure('source piechart')
    explode = (0.1, 0, 0)
    if isOpenmut:
        labels='Exclusively ONT transcriptome','Exclusively reference (gencode)', 'Both'
    else:
        labels='Novel','Annotated','Combination'
    plt.pie([len(ont_only),len(ref_only),len(both)],autopct='%1.1f%%', explode=explode,colors=['#de2d26','#3182bd','#756bb1'])
    plt.title('Peptide spectral hits by source',fontsize=35)
    plt.legend(labels,loc=8)
    plt.savefig(figname)
    plt.clf()
    return("saved to sources_spectral_hits")

def plot_chromosomal_dist(distr_classic,distr_openmut):
    # sns.set(rc={'figure.figsize':(11.7,8.27)})
    # sns.set_style(style='white')
    # combi={'Classical proteogenomics':distr_classic,'Open mutation':distr_openmut}
    plt.figure('chromosomal distribution')
    # combi_df=pd.DataFrame(combi).stack().reset
    # index=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","unknown"]
    chist=pd.DataFrame.from_dict(distr_classic,orient='index').sort_index()
    chist_openmut=pd.DataFrame.from_dict(distr_openmut,orient='index').sort_index()
    combi=pd.concat([chist,chist_openmut],axis=1)
    combi.columns=['Classical proteogenomics','Open mutation']
    combi.plot(kind='bar',legend=False,title="Chromosomal distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Chromosomes")
    plt.savefig('chromosomal_distribution.png')
    plt.clf()
    return("plotted chromosomal distribution")

def plot_coverage_plots(cpdt_pep,fullseqs,fignamehorizontal,fignamevertical):
    '''this function will plot the graphs that correspond to the coverage of the proteome
    - vertical coverage
    - horizontal coverage
    - chromosome distribution
    '''
    high_cov_hor,cov_vert,perc_cov_dist=coverage_measure(cpdt_pep,fullseqs)
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    #horizontal coverage
    plt.figure('horizontal coverage')
    plt.hist(perc_cov_dist,bins=300)
    plt.xlim(1,100)
    plt.ylim(0,1000)
    plt.xlabel("% Coverage")
    plt.ylabel("# Proteins")
    plt.savefig(fignamehorizontal)
    plt.clf()
    #vertical coverage
    plt.figure('vertical coverage')
    plt.hist(cov_vert,bins=1000)
    plt.xlim(1,50)
    plt.xlabel("# Proteins")
    plt.ylabel("# Peptides")
    plt.savefig(fignamevertical)
    plt.clf()
    return("Plotted coverage")

def calc_mut_abundances(mutant_cpdtpep,cpdtpep):
    prot_abundance=[]
    #nr_mutant=[]
    mut_proteins_detected=set()
    num_peptides=0
    num_occurences=0
    for prot,peps in mutant_cpdtpep.items():
        if '_h' in prot:
            stem=prot.split('_h')[0]
        else:
            stem=prot
        if stem in cpdtpep:
            sum_mut=0 #total number of detected mutant peptides
            sum_nonmut=0
            for pep,ct in peps.items():
                sum_mut+=ct
                if ct>0:
                    mut_proteins_detected.add(prot)
                    num_occurences+=ct
                    num_peptides+=1
            if sum_mut>0: #only record the proteins with at least 1 detected mutation peptide
                #nr_mutant.append(sum_mut)
                for pepc,ctc in cpdtpep[stem].items():
                    sum_nonmut+=ctc
                prot_abundance.append((sum_nonmut,sum_mut))
                #prot_abundance.append(sum_mut+sum_nonmut)
    print("Total of "+str(num_occurences)+" occurances of "+str(num_peptides)+" peptides from "+str(len(mut_proteins_detected))+" proteins were detected")
    return(prot_abundance)

def plot_mut(mutant_cpdtpep,cpdtpep,figname):
    '''plot protein abundance vs number of detected mutant peptides'''
    prot_abundance=calc_mut_abundances(mutant_cpdtpep,cpdtpep)
    #make plot
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('mutant peptides')
    plt.scatter(*zip(*prot_abundance))
    plt.xlabel('Protein abundance')
    plt.ylabel('Number mutant peptides detected')
    plt.title('Mutant peptide abundance vs total non-mutated protein abundance')
    plt.savefig(figname)
    plt.clf()
    return('done')

def plot_final_venns(mut_peptide_dict_classic,mut_peptide_dict_openmut,mut_cpdt_theoretical,mutprotset):
    allmuts_classic=Counter()
    allmuts_openmut=Counter()
    #go through each dictionary and concatenate counters so one large counter object with peptides
    for prot,mutct in mut_peptide_dict_classic.items(): #no open mutation
        allmuts_classic+=mutct
    for prott,muti in mut_peptide_dict_openmut.items():
        allmuts_openmut+=muti
    #create diagrams
    plt.figure('venn mutant peptides')
    vda=venn2_unweighted([allmuts_classic,allmuts_openmut],('Proteogenomics approach','Open mutation search')) #venn for the overlap in detected peptides
    plt.title("Unique observed variant peptides",fontsize=26)
    for text in vda.set_labels:
        text.set_fontsize(26)
    for text in vda.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_peps.png')
    plt.clf()
    plt.figure('venn mutant proteins')
    vdb=venn2_unweighted([mut_peptide_dict_classic.keys(),mut_peptide_dict_openmut.keys()],("Proteogenomics approach","Open mutation search")) #venn for the overlap in detected proteins
    plt.title("Unique proteins associated with observed variant peptides",fontsize=26)
    for text in vdb.set_labels:
        text.set_fontsize(26)
    for text in vdb.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_prots.png')
    plt.clf()
    plt.figure('venn proteins all')
    vdb=venn3([mut_peptide_dict_classic.keys(),mut_cpdt_theoretical.keys(),list(mutprotset)],("Proteogenomics approach","All theoretical","All predicted open mutation")) #venn for the overlap in detected proteins
    plt.title("Unique proteins associated with mutations",fontsize=26)
    for text in vdb.set_labels:
        text.set_fontsize(26)
    for text in vdb.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_all_detected_mut_prots.png')
    plt.clf()
    return('plotted final venns')

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

def detect_mut_peptides(pep,ids,cpdt_pep):
    for isi in ids:
        i=get_id(isi)
        found=False
        possibilities=[i,i+'_h0',i+'_h1']
        for poss in possibilities:
            if poss in cpdt_pep.keys():
                for p in cpdt_pep[poss]:
                    if pep in p:
                        return(poss)
                # if pep in cpdt_pep[poss]:
                #     return(poss)
    return('')

def add_to_observed_mutdict(mut_prot,pep,olddict):
    newdict=olddict
    if mut_prot not in newdict:
        newdict[mut_prot]=Counter()
    newdict[mut_prot][pep]+=1
    return(newdict)
    
def combidict_analysis(combidict,chromdict,cpdt_pep,full_seqs,mut_cpdt_theoretical,isOpenmut):
    proteins_covered=Counter() #proteins detected
    mutated=set() #all proteins that were detected to have a mutation by ionbot. how does compare to the proteins that actually do have mutation?
    ref_only=set() #scan ids in the reference set
    ont_only=set() #scan ids in the ont set
    both=set() #scan ids that matched to both ref and ont proteins
    hits_missed=0
    hit_mut=0
    hits_missed_mut=0
    chrom_dist=Counter() # 1 chromosome location per scan id
    mut_cpdt_observed={}
    mutprotset=set()
    for row in tqdm(combidict.iterrows()):
        scanid=row[1][0]
        mod=str(row[1][7])
        aamod=re.findall('[A-Z]->[A-Z]',mod)
        pep=row[1][3]
        if '||' in row[1][9]:
            ids=row[1][9].split('||')
        else:
            ids=[row[1][9]]
        proteins_covered=detected_proteins(ids,proteins_covered) #what proteins from the proteome are covered and in what amounts
        cpdt_pep,notfound=fill_cpdt(pep,mod,ids,cpdt_pep) #what peptides are detected, how many, and what proteins they come from
        hits_missed+=notfound
        chrom_origin=find_chrom(ids,chromdict)
        chrom_dist[chrom_origin]+=1 #which chromosome does the peptide belong to
        ref_only,ont_only,both=bin_hits_by_source(scanid,ids,ref_only,ont_only,both, isOpenmut)
        if len(aamod)>0: #if ib detects a mutated peptide (for open mutation search only)
            mutprotset=mutprotset.union(set(ids))
            hit_mut+=1
            #mut_cpdt_pep,notfound_mut=fill_cpdt()
            mut_prot=detect_mut_peptides(pep,ids,mut_cpdt_theoretical)
            #add mutant peptide to observed
            if mut_prot!='':
                mut_cpdt_observed=add_to_observed_mutdict(mut_prot,pep,mut_cpdt_observed)
            #hits_missed_mut+=notfound_mut
            for i in ids: 
                mutated.add(get_id(i))
            #mutdict_id,pep
        elif not isOpenmut: #check if mutant peptide if not open mutation settings
            mutcand=detect_mut_peptides(pep,ids,mut_cpdt_theoretical)
            if mutcand!='': 
                if pep in mut_cpdt_theoretical[mutcand]:
                    mut_cpdt_observed=add_to_observed_mutdict(mutcand,pep,mut_cpdt_observed)
    #create the figures
    print("number of hits with detected mutation = " +str(hit_mut)+ " matched to "+str(len(mutated))+ " proteins.")
    print("number of hits that were not counted because they were not predicted by in silico digest: "+str(hits_missed))
    print("number of mutant peptides not matched to predicted mutant peptides = " +str(hits_missed_mut))
    #create checkpoint- save the results from above so that whole analysis does not need to be repeated to re-create the graphs
    if isOpenmut:
        # with open('checkpoint.openmut.json','w'):
        #     json.dump()
        plot_mut(mut_cpdt_observed,cpdt_pep,"mutant_abundance_varfree.png")
        plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varfree.png","vertical_coverage_varfree.png")
        plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varfree.png",isOpenmut)
        # plot_chromosomal_dist(chrom_dist,"chromosomal_distribution_varfree.png")
    else:
        plot_mut(mut_cpdt_observed,cpdt_pep,"mutant_abundance_varcont.png")
        plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varcont.png","vertical_coverage_varcont.png")
        plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varcont.png",isOpenmut)
        # plot_chromosomal_dist(chrom_dist,"chromosomal_distribution_varcont.png")
    if isOpenmut:
        return(mut_cpdt_observed,mutprotset,chrom_dist)
    return(mut_cpdt_observed,chrom_dist)


def create_chromosome_reference(gfffile,bedfile):
    chromdict_ref=import_gff(gfffile,False) #import gff3 file annotations from gencode
    chromdict_ont=import_gff(bedfile,True) #import bed file annotations from ont (converted from psl)
    chromdict={**chromdict_ont,**chromdict_ref} #combine the 2 dictionaries
    return(chromdict)

def main(directory_ontonly, directory_refonly, directory_combination, directory_combination_including_variants, cpdtfile,cpdtfile_mut,gfffile,bedfile):
    '''this is the main function that will iterate over the giant pandas df and perform all analyses and make all figures
    this is written in a way that iteration should only be done once.
    '''
    #imports
    print("importing data")
    ibdf_ontonly=concatenate_csvs(directory_ontonly)
    ibdf_refonly=concatenate_csvs(directory_refonly)
    ibdf_combi=concatenate_csvs(directory_combination)
    ibdf_combi_pg=concatenate_csvs(directory_combination_including_variants)

    #qc function
    plot_scores(ibdf_ontonly.dropna(),ibdf_refonly.dropna(),ibdf_combi.dropna())
    plot_scores_pg(ibdf_combi.dropna(),ibdf_combi_pg.dropna())

    #filter badly scoring hits
    ibdf_ontonly = chunk_preprocessing(ibdf_ontonly)
    ibdf_refonly = chunk_preprocessing(ibdf_refonly)
    ibdf_combi = chunk_preprocessing(ibdf_combi)
    ibdf_combi_pg= chunk_preprocessing(ibdf_combi_pg)

    #import insilico digest info
    cpdt_pep,full_seqs=import_cpdt(cpdtfile) #import cpdt will all peptides (cat gencode and flair beforehand). full seqs for calculating horizontal coverage
    mut_cpdt_theoretical=import_cpdt_simple(cpdtfile_mut) #import the cpdt file with all snv peptides
    chromdict=create_chromosome_reference(gfffile,bedfile) #import information about the chromosome of origin (QC)
    
    #iterate to fill the data structures
    print("Analyzing data...")
    mut_observed_openmut,mutprotset,chromdist_openmut=combidict_analysis(ibdf_combi,chromdict,cpdt_pep,full_seqs,mut_cpdt_theoretical,True)
    plt.clf()
    mut_observed_classic,chromdist_classic=combidict_analysis(ibdf_combi_pg,chromdict,cpdt_pep,full_seqs,mut_cpdt_theoretical,False)
    plt.clf()
    plot_chromosomal_dist(chromdist_classic,chromdist_openmut)
    plot_final_venns(mut_observed_classic,mut_observed_openmut,mut_cpdt_theoretical,mutprotset)
    return("Finished")
    
    
main(*sys.argv[1:]) #args: folder containing ONT-variant-free ionbot search, folder containing ref-variant-free ionbot search, folder containing combi-variant-free ionbot search results, combi-variant-free peptides, combi-variant-containing peptides
    
    
    
    
    