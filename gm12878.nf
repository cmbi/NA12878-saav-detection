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
 * PART 1: Create variant-free dictionary
 *
 * Process 1A: Identify coding sequences by running them through ANGEL and SQANTI and then identifying the overlap
 */
process 'identify_coding' { 

    input:
        file 
}

/*
* Process 1B: combine the novel sequences verified by step 1A with the reference proteome (pc_translations)
*/



/**********
 * PART 2: Create variant-containing sequences
 *
 * Process 2A: Filter gff for CDS sequences only
 */

process 'filter_gff' { 

  input: 
      file gff_file
  output:
      file 'dbs/db_fa' into search_db_path
      file 'dbs/db_reversed_fa' into search_db_path_decoy
  
  script:
  """
  grep 'CDS' $gff_file > $bed_file
  """
}

/*
* Process 2B: Get corresponding genomic fasta sequences with bedtools getfasta (both ref and sample)
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
  python3 ../../../gm12878/customized_sequence_creation/translate_annotated.py --fasta cds.variantcontaining.fa --report gencode.v29.vartracking.cds.txt --pfasta gencode.v29.variant_containing_proteins.cds_based.fa

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

