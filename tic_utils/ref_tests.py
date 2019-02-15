from __future__ import print_function
from Bio import SeqIO
import os
import matplotlib.pyplot as plt

def ref_search(fasta,name,seq,n):
	#---Search sequences for characteristic motifs and collect statistics
	has_Nterm = has_N_term(seq)
	has_Intern = has_motif(seq)
	nref = 0
	fasta = fasta.split('.')[0]+'_'
	if n == 1:
		open(fasta+'has_motifs.fasta','w').close()
		open(fasta+'partial_motifs.fasta','w').close()
		open(fasta+'no_motifs.fasta','w').close()
	elif has_Nterm and has_Intern:
		nref = 1
		with open(fasta+'has_motifs.fasta','a') as f:
			f.write(">"+name+"\n"+seq+"\n")
	elif not has_Nterm and not has_Intern:
		with open(fasta+'no_motifs.fasta','a') as f:
			f.write(">"+name+"\n"+seq+"\n")
	else:
		with open(fasta+'partial_motifs.fasta','a') as f:
			f.write(">"+name+"\n"+seq+"\n")
	return nref

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
