import sys
import os
import pandas as pd
from Bio import SeqIO

#create concatenated fasta for percolator
command = "cp %s fasta.tmp"%sys.argv[2]
os.system(command)
with open("fasta.tmp","a") as f:
	for record in SeqIO.parse(sys.argv[2], "fasta"):
	    f.write(">Random_"+record.id+"\n")
            f.write(str(record.seq[::-1])+"\n")

def getLabel(x):
	if x=="T": return 1
	return -1

d = pd.read_csv(sys.argv[1])

tmp = pd.DataFrame()
tmp["Spec_Id"] = d["scan_id"]
tmp["Label"] = d["DB"].apply(getLabel)
tmp["ScanNr"] = [i for i in range(len(d))]
tmp["charge"] = d["charge"]
tmp["F0"] = d["percolator_psm_score"]
tmp["Peptide"] = "-."+ d["matched_peptide"].astype(str)+".-"
tmp["Proteins"] = d["proteins"]
tmp.to_csv("r.pin",sep="\t",index=False)

command = "sed 's/||/\t/g' r.pin > s.pin"
os.system(command)
command = "percolator -i 10 -U -N 200000 -P Random_ -f %s -g -l result.proteins.tsv s.pin > result.out" % ("fasta.tmp")
os.system(command)
