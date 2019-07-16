#!/usr/bin/env python3

import matplotlib, re, os,sys, itertools
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
            temp["scan_id"]=temp["scan_id"].astype(str)+'_'+csvname #make scan ids unique again when concatenating all files
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
    stranddict={}
    with open(gfffile) as handle:
        for line in handle:
            info=line.split('\t')
            if len(info)>5:
                if not isBed:
                    if 'transcript_type=protein_coding' in line:
                        tid=line.split('transcript_id=')[1]
                        tid=tid.split(';')[0]
                        chromdict[tid]=info[0]
                        stranddict[tid]=info[6]
                else:
                    chromdict[info[3]]=info[0]
                    if info[5]!='+' and info[5]!='-':
                        stranddict[info[3]]='unknown'
                    else:
                        stranddict[info[3]]=info[5]
    return(chromdict,stranddict)

def find_chrom(prots,chromdict):
    for p in prots:
        p=get_id(p)
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
        p=get_id(p)
        if '_h' in p:
            p=p.split('_h')[0]
        if p in stranddict:
            return(stranddict[p])
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
        if len(peps.keys())>0:
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
            vert=float(count_pep/len(peps.keys()))
            if vert>0:
                vert_cov.append(vert)
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

def plot_support(prot_evidence,unamb_prot_evidence,figname):
    '''look into the support for proteins
    how many proteins have 1,2,3... peptides supporting their existence? unambiguously?
    1 counter that counts all peptides, another counter that counts only unique peptides'''
    plt.figure('support')
    recountpev=counter_translator(prot_evidence)
    recountunamb=counter_translator(unamb_prot_evidence)
    new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,'20+']
    pev= pd.DataFrame.from_dict(recountpev,orient='index')
    punamb= pd.DataFrame.from_dict(recountunamb,orient='index')
    combi=pd.concat([pev,punamb],axis=1,sort=True)
    combi=combi.reindex(new_index)
    combi.fillna(0,inplace=True)
    combi.columns=['Any peptide evidence','Unambiguous assignments']
    combi.plot(kind='bar',legend=False,title="Peptide support for proteins")
    plt.ylabel("# proteins")
    plt.ylim(0,20000)
    plt.xlabel("# peptides")
    plt.legend(loc='upper right')
    plt.savefig(figname)
    plt.clf()
    return('Plotted peptide support')

def counter_translator(counterobj):
    ct_prots=Counter()
    for elem,ct in counterobj.items():
        if ct>19:
            ct_prots['20+']+=1
        else:
            ct_prots[ct]+=1
    return(ct_prots)

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
    sns.distplot(ibdf_refonly['ionbot_psm_score'], hist=False, label='Human reference only',axlabel='Ionbot score')
    sns.distplot(ibdf_ontonly['ionbot_psm_score'], hist=False, label='Transcriptome translation only',axlabel='Ionbot score')
    sns.distplot(ibdf_combi['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation',axlabel='Ionbot score')
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
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='T']['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation (no variants)',axlabel='Ionbot score')
    sns.distplot(ibdf_combi_pg[ibdf_combi_pg['DB']=='T']['ionbot_psm_score'], hist=False, label='Reference + transcriptome translation (with variants)',axlabel='Ionbot score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig("qc_pearsonr_pgvsopenmut.png")
    plt.clf()
    return("Scores plot made")

def plot_scores_decoy(ibdf_combi,figname):
    '''look at the quality of the matches per dictionary before the dataset has been filtered'''
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure("Pearson R distribution 2")
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='T']['ionbot_psm_score'], label='Target',axlabel='Ionbot score')
    sns.distplot(ibdf_combi[ibdf_combi['DB']=='D']['ionbot_psm_score'], label='Decoy',axlabel='Ionbot score')
    plt.legend()
    plt.title('Correlation between theoretical and observed spectra of matched peptide')
    plt.savefig(figname)
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
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M','unknown']
    plt.figure('chromosomal distribution')
    chist=pd.DataFrame.from_dict(distr_classic,orient='index')#.sort_index()
    chist_openmut=pd.DataFrame.from_dict(distr_openmut,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_openmut],axis=1)
    combi=combi.reindex(new_index)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi.plot(kind='bar',legend=False,title="Chromosomal distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Chromosome")
    plt.legend(loc='upper right')
    plt.savefig('chromosomal_distribution.png')
    plt.clf()
    return("plotted chromosomal distribution")

def plot_strand_dist(distr_classic,distr_openmut):
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('strand distribution')
    chist=pd.DataFrame.from_dict(distr_classic,orient='index')#.sort_index()
    chist_openmut=pd.DataFrame.from_dict(distr_openmut,orient='index')#.sort_index()
    combi=pd.concat([chist,chist_openmut],axis=1)
    combi.columns=['Combi variant-containing','Combi variant-free']
    combi.plot(kind='bar',legend=False,title="Strand distribution of peptide hits")
    plt.ylabel("# Peptides")
    plt.xlabel("Strand")
    plt.legend(loc='upper right')
    plt.savefig('strand_distribution.png')
    plt.clf()
    return("plotted strand distribution")

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
    plt.title("Horizontal coverage")
    plt.xlabel("% Coverage")
    plt.ylabel("# Proteins")
    plt.savefig(fignamehorizontal)
    plt.clf()
    #vertical coverage
    plt.figure('vertical coverage')
    plt.hist(cov_vert,bins=80)
    plt.xlim(1,40)
    plt.ylim(0,4000)
    plt.title("Vertical coverage")
    plt.xlabel("Peptide count (protein size normalized)")
    plt.ylabel("Density")
    plt.savefig(fignamevertical)
    plt.clf()
    return("Plotted coverage")

def calc_nsaf_standard(cpdt_pep,fullseqs):
    '''after cpdt_pep has been filled, find the sum nsaf in order to standardize the abundance scores in calc_mut_abundances'''
    nsaf=0
    for prot,peps in cpdt_pep.items():
        count_peps=0
        for p,ct in peps.items():
            count_peps+=ct
        length_prot=len(fullseqs[prot])
        nsaf+=float(count_peps/length_prot)
    return(nsaf)

def calc_nsaf_protein(pepdict_singleprot,lenprot,sumnsaf):
    sum_nonmut=0
    for normpep,normct in pepdict_singleprot.items():
        sum_nonmut+=normct
    nsaf=float(float(sum_nonmut/lenprot)/sumnsaf)
    return(nsaf)

def calc_mut_abundances(mutant_cpdtpep,cpdtpep,fullseqs):
    mut_pep_abundance=[]
    nonmut_pep_abundance=[]
    #nr_mutant=[]
    mut_proteins_detected=set()
    num_peptides=0
    num_occurences=0
    sumnsaf=calc_nsaf_standard(cpdtpep,fullseqs)
    for prot,peps in mutant_cpdtpep.items():
        if '_h' in prot:
            stem=prot.split('_h')[0]
        else:
            stem=prot
        if stem in cpdtpep:
            nsafnonmut=calc_nsaf_protein(cpdtpep[stem],len(fullseqs[stem]),sumnsaf)
            sum_mut=0 #total number of detected mutant peptides
            sum_nonmut=0
            for pep,ct in peps.items():
                sum_mut+=ct
                if ct>0:
                    mut_proteins_detected.add(prot)
                    num_occurences+=ct
                    num_peptides+=1
                    mut_pep_abundance.append((nsafnonmut,ct))
            #calculate count of non-mutant peptides
            # for normpep,normct in cpdtpep[stem].items():
            #     sum_nonmut+=normct
            # lennonmut=len(fullseqs[stem])
            for normpep,normct in cpdtpep[stem].items():
                # sum_nonmut+=normct
                nonmut_pep_abundance.append((nsafnonmut,normct))
            # nsafnonmut=float(float(sum_nonmut/lennonmut)/sumnsaf)
            # nonmut_pep_abundance.append((nsafnonmut,sum_nonmut))
            # if sum_mut>0: #only record the proteins with at least 1 detected variant peptide
            #     #nr_mutant.append(sum_mut)
            #     mut_pep_abundance.append((nsafnonmut,sum_mut))
            #     #mut_pep_abundance.append(sum_mut+sum_nonmut)
    print("Total of "+str(num_occurences)+" occurances of "+str(num_peptides)+" peptides from "+str(len(mut_proteins_detected))+" proteins were detected")
    return(mut_pep_abundance,nonmut_pep_abundance)

def calculate_correlation(mut_pep_abundance,nonmut_pep_abundance):
    dt=np.dtype('float,int')
    variant = np.array(mut_pep_abundance,dtype=dt)
    nonvariant=np.array(nonmut_pep_abundance,dtype=dt)
    cor_var= np.corrcoef(variant['f0'],variant['f1'])
    cor_nonvar= np.corrcoef(nonvariant['f0'],nonvariant['f1'])
    # cor_var_nonvar=np.corrcoef(variant['f0'],nonvariant['f0'])
    print("correlation between variant peptides abundance and total peptide abundance: "+str(cor_var[1][0]))
    print("correlation between non-variant peptide abundance and total peptide abundance: "+str(cor_nonvar[1][0]))
    return(0)

def plot_mut(mutant_cpdtpep,cpdtpep,fullseqs,figname):
    '''plot protein abundance vs number of detected mutant peptides'''
    mut_pep_abundance,nonmut_pep_abundance=calc_mut_abundances(mutant_cpdtpep,cpdtpep,fullseqs)
    calculate_correlation(mut_pep_abundance,nonmut_pep_abundance)
    #make plot
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style(style='white')
    plt.figure('mutant peptides')
    plt.scatter(*zip(*mut_pep_abundance),c='r',label='Mutant peptide',alpha=0.5)
    plt.scatter(*zip(*nonmut_pep_abundance),c='b',label='Normal peptide',alpha=0.5)
    plt.xlabel('Protein abundance (NSAF normalized)')
    plt.xlim(0,0.001)
    plt.ylabel('Number mutant peptides detected')
    plt.title('Peptide abundance vs total protein abundance')
    plt.legend(loc='upper right')
    plt.savefig(figname)
    plt.clf()
    return('done')

def discrepancy_check(mut_peptide_dict_classic,mut_peptide_dict_openmut,ibdf_combi,ibdf_combi_pg):
    '''check out the differences in identifications between the 2 combination dictionaries
    why doesn't ionbot catch everything? look at the ones that it does not catch but the variant-containing dictionary does
    plot lengths of the missed/caught peptides (longer than average?)
    plot unexpected modifications of the missed peptides (more unexpected modifications than average?)'''
    allmuts_classic=Counter()
    allmuts_openmut=Counter()
    for prot,mutct in mut_peptide_dict_classic.items(): #no open mutation
        allmuts_classic+=mutct
    for prott,muti in mut_peptide_dict_openmut.items():
        allmuts_openmut+=muti
    discrepancy=set(allmuts_classic).difference(set(allmuts_openmut))
    discrepancy_ct_pg=allmuts_classic-allmuts_openmut
    discrepancy_ct_om=allmuts_openmut-allmuts_classic
    #plot lengths of the peptides caught by one method but not the other
    plot_peplengths(discrepancy_ct_pg,discrepancy_ct_om)
    #get the scan ids corresponding to all the variant peptides that are in the variant containing search results but not variant free
    scanids=ibdf_combi_pg.loc[ibdf_combi_pg["matched_peptide"].isin(discrepancy),"scan_id"].tolist()
    #return dataframe containing only rows from other result dictionary corresponding to the list of scan ids just obtained
    discr_df=ibdf_combi.loc[ibdf_combi["scan_id"].isin(scanids)]
    plot_unexpected_mods(discr_df["unexpected_modification"].tolist())
    return(0)

def plot_unexpected_mods(list_mods):
    mod_ct=categorize_mods(list_mods)
    plt.figure('discrepant peptide lengths')
    chist_pg=pd.DataFrame.from_dict(list_mods.most_common(20),orient='index')
    chist_pg.plot(kind='bar',legend=False,title="Unexpected modifications found instead of SAAVs from variant peptides (variant-free search)")
    plt.ylabel("Density")
    plt.xlabel("Length peptide")
    plt.legend(loc='upper right')
    plt.savefig('discrepant_peptide_length.png')
    plt.clf()
    return(0)

def categorize_mods(list_mods):
    mod_ct=Counter()
    for mod in list_mods:
        if len(re.findall('[A-Z]->[A-Z]',str(mod)))>0:
            mod_ct['SAAV']+=1
        elif mod=='nan':
            mod_ct["none"]+=1
        else:
            s_mod=re.split('\[[a-z]\]',mod)[0]
            mod_ct[s_mod]+=1
    return(mod_ct)

def plot_peplengths(peptide_counter_pg,peptide_counter_om):
    lenct_pg=gather_counts(peptide_counter_pg)
    lenct_om=gather_counts(peptide_counter_om)
    plt.figure('discrepant peptide lengths')
    # new_index= [1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M','unknown']
    chist_pg=pd.DataFrame.from_dict(peptide_counter_pg,orient='index').sort_index()
    chist_om=pd.DataFrame.from_dict(peptide_counter_om,orient='index').sort_index()
    combi=pd.concat([chist_pg,chist_om],axis=1)
    combi.columns=['Combi variant-containing only','Combi variant-free only']
    combi.plot(kind='bar',legend=False,title="Length of discrepant peptides (found by one method and not the other)")
    # combi=combi.reindex(new_index)
    plt.ylabel("Density")
    plt.xlabel("Length peptide")
    plt.legend(loc='upper right')
    plt.savefig('discrepant_peptide_length.png')
    plt.clf()
    return(0)

def gather_counts(peptide_counter):
    length_counter=Counter()
    for pep in peptide_counter:
        length_counter[len(pep)]+=peptide_counter[pep]
    return(length_counter)

def plot_final_venns(mut_peptide_dict_classic,mut_peptide_dict_openmut,mut_cpdt_theoretical,mutprotset):
    allmuts_classic=Counter()
    allmuts_openmut=Counter()
    #go through each dictionary and concatenate counters so one large counter object with peptides
    for prot,mutct in mut_peptide_dict_classic.items(): #no open mutation
        allmuts_classic+=mutct
    for prott,muti in mut_peptide_dict_openmut.items():
        allmuts_openmut+=muti
    # print(set(allmuts_classic).difference(set(allmuts_openmut))) #TESTING PURPOSES
    #create diagrams
    plt.figure('venn mutant peptides')
    vda=venn2_unweighted([allmuts_classic,allmuts_openmut],('Combi variant-containing','Combi variant-free')) #venn for the overlap in detected peptides
    plt.title("Unique observed variant peptides",fontsize=26)
    for text in vda.set_labels:
        text.set_fontsize(26)
    for text in vda.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_peps.png')
    plt.clf()
    plt.figure('venn mutant proteins')
    vdb=venn2_unweighted([mut_peptide_dict_classic.keys(),mut_peptide_dict_openmut.keys()],("Combi variant-containing","Combi variant-free")) #venn for the overlap in detected proteins
    plt.title("Unique proteins associated with observed variant peptides",fontsize=26)
    for text in vdb.set_labels:
        text.set_fontsize(26)
    for text in vdb.subset_labels:
        text.set_fontsize(20)
    plt.savefig('overlap_detected_mut_prots.png')
    plt.clf()
    plt.figure('venn proteins all')
    vdb=venn3([set(mut_peptide_dict_classic.keys()),set(mut_cpdt_theoretical.keys()),mutprotset],("Combi variant-containing","All theoretical","Combi variant-free")) #venn for the overlap in detected proteins
    plt.title("Unique proteins associated with variant peptides",fontsize=26)
    # for text in vdb.set_labels:
    #     text.set_fontsize(26)
    # for text in vdb.subset_labels:
    #     text.set_fontsize(20)
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

def detect_mut_peptides(pep,ids,cpdt_pep,isOpenmut):
    for isi in ids:
        i=get_id(isi)
        found=False
        if isOpenmut:
            possibilities=[i,i+'_h0',i+'_h1']
        else:
            possibilities=[i]
        for poss in possibilities:
            if poss in cpdt_pep.keys():
                for p in cpdt_pep[poss]: #this allows for some flexibility in matches, allows for non-canonical
                    if pep in p or equivalent_check(pep,p):
                        return(poss)
    return('')

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

def add_to_observed_mutdict(mut_prot,pep,olddict):
    newdict=olddict
    if mut_prot not in newdict:
        newdict[mut_prot]=Counter()
    newdict[mut_prot][pep]+=1
    return(newdict)
    
def combidict_analysis(combidict,chromdict,stranddict,cpdt_pep,full_seqs,mut_cpdt_theoretical,isOpenmut):
    proteins_covered=Counter() #proteins detected
    mutated=set() #all proteins that were detected to have a variant by ionbot. how does compare to the proteins that actually do have variant?
    ref_only=set() #scan ids in the reference set
    ont_only=set() #scan ids in the ont set
    both=set() #scan ids that matched to both ref and ont proteins
    hits_missed=0
    hit_mut=0
    hits_missed_mut=0
    chrom_dist=Counter() # 1 chromosome location per scan id
    strand_dist=Counter()
    mut_cpdt_observed={}
    protein_support=Counter()
    unamb_protsupport=Counter()
    for row in tqdm(combidict.iterrows()):
        scanid=row[1][0]
        mod=str(row[1][7])
        aamod=re.findall('[A-Z]->[A-Z]',mod)
        pep=row[1][3]
        if '||' in row[1][9]:
            ids=row[1][9].split('||')
            for i in ids:
                protein_support[get_id(i)]+=1
        else: #unambiguous assignment!
            ids=[row[1][9]]
            unamb_protsupport[get_id(row[1][9])]+=1
        proteins_covered=detected_proteins(ids,proteins_covered) #what proteins from the proteome are covered and in what amounts
        cpdt_pep,notfound=fill_cpdt(pep,mod,ids,cpdt_pep) #what peptides are detected, how many, and what proteins they come from
        hits_missed+=notfound
        chrom_origin=find_chrom(ids,chromdict)
        strand_origin=find_strand(ids,stranddict)
        chrom_dist[chrom_origin]+=1 #which chromosome does the peptide belong to
        strand_dist[strand_origin]+=1 #which strand does the peptide belong to
        ref_only,ont_only,both=bin_hits_by_source(scanid,ids,ref_only,ont_only,both, isOpenmut)
        if isOpenmut and len(aamod)>0: #if ib detects a mutated peptide (for open variant search only)
            hit_mut+=1
            #mut_cpdt_pep,notfound_mut=fill_cpdt()
            mut_prot=detect_mut_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)
            #add mutant peptide to observed
            if mut_prot!='':
                mut_cpdt_observed=add_to_observed_mutdict(mut_prot,pep,mut_cpdt_observed)
            #hits_missed_mut+=notfound_mut
            for i in ids: 
                mutated.add(get_id(i))
            #mutdict_id,pep
        # elif isOpenmut and detect_mut_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)!='': ##very strange scenario here!!##
        #     print(scanid)
        elif not isOpenmut: #check if mutant peptide if not open mutation settings
            mutcand=detect_mut_peptides(pep,ids,mut_cpdt_theoretical,isOpenmut)
            if mutcand!='': 
                if pep in mut_cpdt_theoretical[mutcand]:
                    mut_cpdt_observed=add_to_observed_mutdict(mutcand,pep,mut_cpdt_observed)
    #create the figures
    print("number of hits with detected variant = " +str(hit_mut)+ " matched to "+str(len(mutated))+ " proteins.")
    print("number of hits that were not counted because they were not predicted by in silico digest: "+str(hits_missed))
    print("number of mutant peptides not matched to predicted mutant peptides = " +str(hits_missed_mut))
    #create checkpoint- save the results from above so that whole analysis does not need to be repeated to re-create the graphs
    if isOpenmut:
        # with open('checkpoint.openmut.json','w'):
        #     json.dump()
        plot_mut(mut_cpdt_observed,cpdt_pep,full_seqs,"mutant_abundance_varfree.png")
        plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varfree.png","vertical_coverage_varfree.png")
        plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varfree.png",isOpenmut)
        plot_support(protein_support,unamb_protsupport,'protein_evidence_varfree.png')
    else:
        plot_mut(mut_cpdt_observed,cpdt_pep,full_seqs,"mutant_abundance_varcont.png")
        plot_coverage_plots(cpdt_pep,full_seqs,"horizontal_coverage_varcont.png","vertical_coverage_varcont.png")
        plot_source_piechart(ref_only,ont_only,both,"sources_spectral_hits_varcont.png",isOpenmut)
        plot_support(protein_support,unamb_protsupport,'protein_evidence_varcont.png')
    if isOpenmut:
        return(mut_cpdt_observed,mutated,chrom_dist,strand_dist)
    return(mut_cpdt_observed,chrom_dist,strand_dist)


def create_chromosome_reference(gfffile,bedfile):
    chromdict_ref,stranddict_ref=import_gff(gfffile,False) #import gff3 file annotations from gencode
    chromdict_ont,stranddict_ont=import_gff(bedfile,True) #import bed file annotations from ont (converted from psl)
    chromdict={**chromdict_ont,**chromdict_ref} #combine the 2 dictionaries
    stranddict={**stranddict_ont,**stranddict_ref} #combine the 2 dictionaries
    return(chromdict,stranddict)

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
    plot_scores_decoy(ibdf_combi.dropna(),"qc_pearsonr_decoy_varfree.png")
    plot_scores_decoy(ibdf_combi_pg.dropna(),"qc_pearsonr_decoy_varcont.png")

    #filter badly scoring hits
    ibdf_ontonly = chunk_preprocessing(ibdf_ontonly)
    ibdf_refonly = chunk_preprocessing(ibdf_refonly)
    ibdf_combi = chunk_preprocessing(ibdf_combi)
    ibdf_combi_pg= chunk_preprocessing(ibdf_combi_pg)

    #import insilico digest info
    cpdt_pep,full_seqs=import_cpdt(cpdtfile) #import cpdt will all peptides (cat gencode and flair beforehand). full seqs for calculating horizontal coverage
    mut_cpdt_theoretical=import_cpdt_simple(cpdtfile_mut) #import the cpdt file with all snv peptides
    chromdict,stranddict=create_chromosome_reference(gfffile,bedfile) #import information about the chromosome of origin (QC)
    
    #iterate to fill the data structures
    print("Analyzing data...")
    mut_observed_openmut,mutprotset,chromdist_openmut,stranddist_openmut=combidict_analysis(ibdf_combi,chromdict,stranddict,cpdt_pep,full_seqs,mut_cpdt_theoretical,True)
    plt.clf()
    mut_observed_classic,chromdist_classic,stranddist_classic=combidict_analysis(ibdf_combi_pg,chromdict,stranddict,cpdt_pep,full_seqs,mut_cpdt_theoretical,False)
    plt.clf()
    discrepancy_check(mut_observed_classic,mut_observed_openmut, ibdf_combi, ibdf_combi_pg)
    plot_chromosomal_dist(chromdist_classic,chromdist_openmut)
    plot_strand_dist(stranddist_classic,stranddist_openmut)
    plot_final_venns(mut_observed_classic,mut_observed_openmut,mut_cpdt_theoretical,mutprotset)
    return("Finished")
    
    
main(*sys.argv[1:]) #args: folder containing ONT-variant-free ionbot search, folder containing ref-variant-free ionbot search, folder containing combi-variant-free ionbot search results, combi-variant-free peptides, combi-variant-containing peptides
    
    
    
    
    