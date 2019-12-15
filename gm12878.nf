#!/usr/bin/env nextflow
    
/*
 * Copyright (c) 
 */
 
  
/* 
 * 
 * Renee Salz
 */


/*
 * Define the default parameters
 */ 

params.db = "$baseDir/data/uniprot_some_crap.fa"
params.missed_cleavages = 1
params.i = "/home/compomics/public3/PJ/proteomics_privacy/OSCC/OSCC-P01/OSCC-P01_raw/*.raw"
params.mgf_store = '/home/compomics/public3/PJ/MGF'
params.ionbot_ptm_file = "$baseDir/data/basic_opt_config.file"
params.ionbot_params = '-e -m'

log.info """\
variant detection pipeline   v 1.0 
======================================
input            : $params.i

search db:
missed_cleavages : $params.missed_cleavages
db               : $params.db

raw file conversion:
mgf_store        : $params.mgf_store

ionbot:
ionbot_ptm_file  : $params.ionbot_ptm_file
ionbot_params    : $params.ionbot_params
"""
variant_gz=
variant_gz_tbi=
gff_file_reference=
psl_file_

df_fasta = file(params.db)
raw_files = file(params.i)
omega_path = params.omega_path
thermo_path = params.thermo_path
mgf_store = params.mgf_store
ionbot_ptm_file = file(params.ionbot_ptm_file)
ionbot_params = params.ionbot_params

/**********
 * PART 0: File preparation
 *
 * Filter the gencode GFF + protein file for completeness (no trunacated seqs)
 */
process 'filter_gencode' { 

    input: 
        file gff3_file
        file pc_translations
    output:
        file 'gencode.v29.complete.cds.gff3' into gff_filtered
        file 'gencode.v29.complete.pc_translations.fa' into pc_translations_filtered

    script:
    """
    #!/usr/bin/env python3

    import collections
    import pandas as pd
    import gffpandas.gffpandas as gffpd
    import phylopandas as ph

    def write_to_fasta(df,outputfile):
        df['id']='>'+df['id']
        df[['id','sequence']].to_csv(outputfile,header=None,index=None,sep='\n')

    gencode=gffpd.read_gff3($gff_file)
    g_proteins=ph.read_fasta($pc_translations)
    g_proteins['id']=g_proteins['id'].str.split('|').apply(lambda x: x[1])
    gencode_df=gencode.df
    gencode_df=gencode_df[gencode_df['type'].isin(['CDS','three_prime_UTR','five_prime_UTR'])]
    gencode_df['id']=gencode_df['attributes'].str.split('Parent=').apply(lambda x: x[1])
    gencode_df['id']=gencode_df['id'].str.split(';').apply(lambda x: x[0])
    check=gencode_df.groupby('id')['type'].agg(lambda x: set(x)).reset_index()
    complete=check.loc[check['type'].str.len()==3,'id']
    write_to_fasta(g_proteins[g_proteins['id'].isin(list(complete))],'gencode.v29.complete.pc_translations.fa')
    gencode_complete=gencode_df[(gencode_df['id'].isin(list(complete))) & (gencode_df['type']=='CDS')]
    gencode.df=gencode_complete.drop('id',axis=1)
    gencode.to_gff3('gencode.v29.complete.cds.gff3')
    """
}
/*
 * Make psl of the ont transcripts into full transcript fasta files 
 * (since FLAIR does not reverse-compliment their minus strand output transcripts)
 */
process 'psl_to_fa' { 

    input:
        file psl
        file reference_genome
    output:
        file "${psl}.fa" into ont_transcript_fasta

    script:
    """
    from pyfaidx import Fasta
    import pandas as pd

    def build_seq(starts,lengths,chromosome,refseq,strand):
        finalseq=''
        for index,num in enumerate(starts):
            subseq=refseq[chromosome][int(num):int(num)+int(lengths[index])]
            if strand=='-':
                finalseq+=str(-subseq)
            else:
                finalseq+=str(subseq)
        return(finalseq)

    reference=Fasta('../../ref_grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa')
    junctions=pd.read_csv('nvrna.flair.isoforms.setA.psl',sep='\t',header=None)
    junctions['tStarts']=junctions[20].str.split(',').apply(lambda x: list(filter(None, x)))
    junctions['tStarts']=junctions.apply(lambda x: x['tStarts'][::-1] if x[8]=='-' else x['tStarts'],axis=1)
    junctions['blocksize']=junctions[18].str.split(',').apply(lambda x: list(filter(None, x)))
    junctions['blocksize']=junctions.apply(lambda x: x['blocksize'][::-1] if x[8]=='-' else x['blocksize'],axis=1)
    junctions['finalseq']=junctions.apply(lambda x: build_seq(x['tStarts'],x['blocksize'],x[13],reference,x[8]),axis=1)
    junctions['header']='>'+junctions[9]
    junctions[['header','finalseq']].to_csv(outputfile,header=None,index=None,sep='\n')
    """
}

/*
 * Run ANGEL and SQANTI to predict coding sequences
 * note: these may need to be done in a docker to have them work properly in the pipeline
 */

process 'run_angel_and_sqanti2' { 

    input:
        file true_transcripts
        file reference_genome
    output:
        folder sqanti_results
        folder angel_results

    script:
    """
    echo 'deb http://security.ubuntu.com/ubuntu xenial-security main' >> /etc/apt/sources.list
    apt-get update
    apt-get install libpng12-0
    python $SQANTI2_directory/sqanti_qc2.py -t $cpu na12878/ont/nvrna.flair.isoforms.setA.true_transcripts.fa ref_grch38/gencode/gencode.v29.annotation.gtf na12878/GRCh38_full_analysis_set_plus_decoy_hla.fa
    python $ANGEL_directory/dumb_predict.py ../../data/custom_seq/ont/nvrna.flair.isoforms.setA.variant_transcripts.nonannotated.fa ../../data/custom_seq/ont/angel_dumb_predict_novel_only/ --cpus $cpu
    """
}

/**********
 * PART 1: Create variant-free dictionary
 *
 * Process 1A: Find overlap of ORF predictors- save to fasta and txt file
 */

process 'find_overlap' { 

    input:
        file 'dumb_reference.final.cds'
        file 'setA_classification.txt'
        file 'setA.faa'
    output:
        file 'setA_classification_angelverified.txt'
        file 'setA_angelverified.faa'

    script:
    """
    #!/usr/bin/env python3

    import re
    import pandas as pd
    import phylopandas as ph

    def write_to_fasta(df,outputfile):
        df['id']='>'+df['id']
        df[['id','sequence']].to_csv(outputfile,header=None,index=None,sep='\n')

    sqantiout=pd.read_csv('setA_classification.txt',sep='\t')
    sqantiout=sqantiout[sqantiout['coding']=='coding']
    angelout=ph.read_fasta('dumb_reference.final.cds')
    sqantiaa=ph.read_fasta('squanti_out/nvrna.flair.isoforms.setA.true_transcripts.renamed_corrected.faa')
    angelout=angelout[~angelout['id'].str.contains('ENST')]
    angelout['id']=angelout['id'].str.split('|').apply(lambda x: x[0])
    sqantiout=sqantiout[sqantiout['isoform'].isin(angelout['id'].unique())]
    sqantiout.to_csv('setA_classification_angelverified.txt',sep='\t',index=False)
    write_to_fasta(sqantiaa[sqantiaa['id'].isin(angelout['id'].unique())].drop_duplicates(subset='sequence',keep='first'),'setA_angel_and_sqanti_verified.faa')
    """
}

/*
* Process 1B: combine the novel sequences verified by step 1A with the reference proteome (pc_translations)
*/

process 'create_final_dict' { 

    input:
        file 'gencode.v29.complete.pc_translations.fa'
        file 'setA_angel_and_sqanti_verified.faa'
    output:
        file 'variant_free_searchdict.fa'

    script:
    """
    cat gencode.v29.complete.pc_translations.fa setA_angel_and_sqanti_verified.faa > variant_free_searchdict.fa
    """
}

/**********
 * PART 2: Create variant-containing sequences
 *
* Process 2A: Get corresponding genomic fasta sequences with bedtools getfasta (both ref and sample)
*/

process 'convert2bed' { 

  input: 
      file gff_file
      file psl_file
  output:
      file 'dbs/db_fa' into search_db_path
      file 'dbs/db_reversed_fa' into search_db_path_decoy
  
  script:
  """
  bedtools getfasta [OPTIONS] -fi $genome -bed $gff_file

  """
}

/*
* Process 2C: incorporate variants from VCF file into the CDS sequences with phase_exons
*/

process 'incorporate_variants_in_cds' { 

  input: 
      file gff_file
      file psl_file
  output:
      file 'dbs/db_fa' into search_db_path
      file 'dbs/db_reversed_fa' into search_db_path_decoy
  
  script:
  """
  python2 ../../../gm12878/customized_sequence_creation/phase_exons.py --cds gencode.v29.cds.fa --vcf ../../variants/NA12878.vcf.gz --out gencode.v29.cds.variant_containing.txt

  """
}

/*
* Process 2D: combine the CDS for every transcript by their rank with combine_exons
*/
 process 'incorporate_variants_in_cds' { 

  input: 
      file gff_file
      file psl_file
  output:
      file 'dbs/db_fa' into search_db_path
      file 'dbs/db_reversed_fa' into search_db_path_decoy
  
  script:
  """
  python3 ../../../gm12878/customized_sequence_creation/combine_exons.py --cpu 40 --exons gencode.v29.cds.variant_containing.txt --jun ../../ref_grch38/gencode/gencode.v29.cds.gff3 --out cds.variantcontaining.fa --report gencode.v29.vartracking.cds.txt

  """
}

/*
* Process 2E: Translate CDSs into protein with translate_annotated
*/
process 'incorporate_variants_in_cds' { 

  input: 
      file gff_file
      file psl_file
  output:
      file 'dbs/db_fa' into search_db_path
      file 'dbs/db_reversed_fa' into search_db_path_decoy
  
  script:
  """
  python3 ../../../gm12878/customized_sequence_creation/translate_annotated.py --fasta cds.variantcontaining.fa --shifts alternative_start_proteins.txt --report gencode.v29.vartracking.cds.txt --pfasta gencode.v29.variant_containing_proteins.cds_based.fa

  """
}

/**********
* PART 3: Make final variant-containing dictionary
*
* Process 3A: replace the sequences in the variant-free dictionary with variant-containing versions where applicable
*/

/**********
* PART 4: Get variant peptides
*
* Process 4A: generate variant peptides with interesting_peptide_finder
*/

