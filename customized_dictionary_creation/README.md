# Customized database creation

All steps to create the dictionaries are recorded within the nextflow file gm12878.nf. Whether this file will run with nextflow has not been tested, but all the components will work if run seperately. An overview of these components is as follows:

Requirements:
	-Python 2.7 (and associated packages in requirements.txt)
	-Bedtools
	-ANGEL (https://github.com/PacificBiosciences/ANGEL)
	-BEDOPS (optional, needed when working with junction files that are not in .bed or .gff format)


1. Filter GENCODE transcripts for only complete sequences
2. Obtain a bed or gff3 file with one exon per line
	- with gff3 files this can be made by filtering the file for only lines that are annotated as "exon"
	- for a BED12 file, each junction in the transcript will have to take it's own line. There is a helper script in misc/bed12_to_bed.py that does this.
	- For PSL, first convert to BED12 with BEDOPS.
3. Run getfasta with bedtools to obtain the exon sequences in fasta format
	- used only with the "-fi", "-fo" and "-bed" functions
4. Run phase_exons.py to incorporate variants
	- arguments: the fasta file created in step 2, gzipped phased vcf file and output file name
5. Run combine_exons.py to create the final transcripts
	- arguments: variant-incorporated exon fasta made in step 3, junction file, 0 or 1 based, and output file name
6. Take the union for combination databases
