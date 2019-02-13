from __future__ import print_function
from Bio import SeqIO
import calcs
import parser
import os
import matplotlib.pyplot as plt

def has_N_term(seq):
	if 'EPM' in seq:
		return True
	else:
		return False

def has_motif(seq):
	if 'PER' in seq and seq.count('MD') > 2:
		return True
	else:
		return False

def is_reflectin(fasta):
	count = 0.0
	for it,s in enumerate(SeqIO.parse(fasta,'fasta')):
		if has_motif(s.seq) or has_N_term(s.seq):
			count+=1
	print(count/it)

if __name__ == "__main__":
	fasta = "all_reflectin/periodic_seqs.fasta"
	is_reflectin(fasta)
