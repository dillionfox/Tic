from __future__ import print_function
from Bio import SeqIO
import Tic
import calcs
import parser
import os
import matplotlib.pyplot as plt

def isref_analysis():
	# reset options list. Fasta will be reset in loop.
	auto_filter = True ; sin_filter = False  ; display = True
	ac_cutoff = -0.2
	dG_scale = 'interface'
	for fil in ['has_motifs.fasta','partial_motifs.fasta','no_motifs.fasta']:
		calcs.MPEX_tools.reset()
		print("working in", fil, "now!")
		try:
			open(fil,'r')
		except IOError:
			print("Run 'is_reflectin' calculation first. Files do not exist") ; exit()
		mpex = Tic.Compute(fasta=fil,dG_scale=dG_scale,standard_plot=True)
		with open(fil.split('.')[0]+'_Tic.fasta','w') as f, open(fil.split('.')[0]+'_nTic.fasta','w') as fn:
			for name,seq,dG,i in parser.parse_fasta(SeqIO,fil,dG_scale):
				calc = mpex.main_calcs(name,seq,dG,i)
				if calc.period != False:
					f.write(">"+calc.name+"\n"+calc.seq+"\n")
				else:
					fn.write(">"+calc.name+"\n"+calc.seq+"\n")
			if auto_filter:
				print(100*float(calcs.MPEX_tools.af)/calcs.MPEX_tools.nchecked, "% are periodic - auto_filter")
			if sin_filter:
				print(100*float(calcs.MPEX_tools.sf)/calcs.MPEX_tools.nchecked, "% are periodic - sin_filter")
			if display:
				if 'DISPLAY' in os.environ:
					plt.show()
				else:
					print("Cannot open display. Saving plot as 'result.png' in working directory")
					plt.savefig('result.png')
	return None

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
