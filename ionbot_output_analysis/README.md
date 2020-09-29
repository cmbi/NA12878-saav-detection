# Analysis of spectral search hits

This directory contains the scripts used to analyze the output of Ionbot search tool

The script interesting_peptide_finder.py should be run to produce the theoretical lists of variant peptides (both target and decoy variant peptide lists must be produced).

The main script to run is ionbot_analysis.py. Run with -h to see all input requirements. 
Note: a retention time prediction file is required, which can be obtainede by running DeepRT.

Example:
python3 ionbot_analysis.py \
--ont gencode_ref/v0.5.0/ont/ \
--ref gencode_ref/v0.5.0/ref/ \
--cvc gencode_ref/v0.5.0/vc \
--cvf gencode_ref/v0.5.0/vf \
--var ionbot_analysis_files/variant_peps_plus_origin.csv \
--decoy ionbot_analysis_files/variant_decoy_peps.csv \
--bed ionbot_analysis_files/nvrna.flair.isoforms.setA.bed \
--gff ionbot_analysis_files/gencode.v29.annotation.gff3 \
--rt ionbot_analysis_files/preds_calibrated.csv \
--crap ionbot_analysis_files/contaminants.fasta