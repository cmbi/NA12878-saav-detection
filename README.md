# gm12878
Project relating to detection of variants on the transcript and protein level of NA12878

This repository contains all the scripts necessary to make customized protein sequences for NA12878 (and theoretically any other cell line) given a phased VCF and a file with annotation of the exon junctions.

Requirements:
	-Python 2.7 (and associated packages in requirements.txt)
	-Bedtools
	-ANGEL (https://github.com/PacificBiosciences/ANGEL)
	-BEDOPS (optional, needed when working with junction files that are not in .bed or .gff format)


The flow is as follows:

1. Obtain a bed or gff3 file with one exon per line
	- with gff3 files this can be made by filtering the file for only lines that are annotated as "exon"
	- for a BED12 file, each junction in the transcript will have to take it's own line. There is a helper script in misc/bed12_to_bed.py that does this.
	- For PSL, first convert to BED12 with BEDOPS.
2. Run getfasta with bedtools to obtain the exon sequences in fasta format
	- used only with the "-fi", "-fo" and "-bed" functions
3. Run phase_exons.py to incorporate variants
	- arguments: the fasta file created in step 2, gzipped phased vcf file and output file name
4. Run combine_exons.py to create the final transcripts
	- arguments: variant-incorporated exon fasta made in step 3, junction file, 0 or 1 based, and output file name
5. (optional) Run the final transcript fasta through misc/prepare_for_angel.py
6. Input the final transcript file to ANGEL (https://github.com/PacificBiosciences/ANGEL) using dumb_predict.py
7. Various analyses for which scripts can be found in the analysis_scripts folder of this repository
