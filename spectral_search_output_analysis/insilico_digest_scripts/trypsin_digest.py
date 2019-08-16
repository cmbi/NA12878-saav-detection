#!/usr/bin/env python3

# import modules used here
import sys
import argparse
from Bio import SeqIO
from pyteomics.parser import cleave
from pyteomics.parser import expasy_rules

def find_str(s, char):
    index = 0

    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index+len(char)] == char:
                    return index

            index += 1

    return -1

# Gather our code in a main() function
def main():
	parser = argparse.ArgumentParser(description='Cleave.', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('input_file', metavar='INPUT', type=str,
										 help='The input file')
	parser.add_argument('-m', dest='missed_cleavages', type=int, default=1, help='Maximum number of missed cleavages')
	parser.add_argument('-out_cpdt',help="cpdt file to write to",required=False)
								


	args = parser.parse_args()
	#print(args.accumulate(args.integers))
	f=open(args.out_cpdt,'w')
	fasta_sequences = SeqIO.parse(open(args.input_file),'fasta')
	# print("{}\t{}\t{}".format('protein', 'sequence', 'start'))

	pep_list = {}
	for fasta in fasta_sequences:
		name, seq = fasta.id, str(fasta.seq)
		pep_list[name] = []
		f.writelines('>'+name+'\n')
		#seq_cut = cleave(seq, expasy_rules['trypsin'], int(args.missed_cleavages))
		seq_cut = cleave(seq, '[KR]', int(args.missed_cleavages))
		#start = 0
		for i,k in enumerate(seq_cut):
			if k in pep_list[name]:
				continue
			pep_list[name].append(k)
			peplen = len(k)

			if len(k) >= 6:
				start = find_str(seq, k)
				# print("{}\t{}\t{}".format(name, k, start))
				f.writelines('\t'+'PEPTIDE '+k+': '+str(start)+'\n')
				if k.startswith('M'):
					f.writelines('\t'+'PEPTIDE '+k+': '+str(start)+'\n')
					# print("{}\t{}\t{}".format(name, k[1:], start+1))

			#start += peplen




# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
	main()
