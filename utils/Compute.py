import main
from main import plt
from table import table
from Bio import SeqIO
import re
import os

acceptable_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

class sequence:
	def __init__(self,name,aaseq):
		self.name = name
		self.aaseq = aaseq

def read_fasta(fasta):
	for s in SeqIO.parse(fasta,'fasta'):
		yield s

def check_seq(seq):
	for aa in seq:
		if aa not in acceptable_aa: 
			return True
	return False

def display(opts):
	if opts['standard_plot']:
		return True
	if opts['fit_to_sin']:
		print 'fit_to_sin'
		return True
	if opts['auto_filter']:
		return True
	return False

def parse_fasta(opts):
	window_size = 19 ; dw = window_size//2 ; 
	for counter,seq in enumerate(read_fasta(opts['fasta'])):
		trim_seq = re.sub('-','', str(seq.seq))
		if check_seq(trim_seq): continue
		RFE = [table.table[s][opts['dG_scale']] for s in trim_seq]
		dG = []
		for ind in range(dw,len(trim_seq)-dw+1):
			s = RFE[ind-dw:ind+dw+1]
			if len(s) != window_size: continue
			dG.append(sum(s))
		yield seq.id,trim_seq,dG,counter

def parse_MPEX(opts):
	nseq = 60
	counter = 0
	dG = []
	rseq = []
	window_size = 19
	dw = window_size//2
	phi_list = []
	for it,line in enumerate(open(opts['MPEX_fil'])):
		if counter >= nseq: break
		try:
			l = line.split('\t') ; first = l[0]
		except:
			continue
		if len(l) == 0: continue
		if l[0] == '"Position"':
			if len(dG) > 0:
				counter += 1
				fullseq = ''.join(seq).upper() 
				yield it,fullseq,dG,it
			dG = [] ; seq = [] ; rseq = [] ; continue
		resname = l[1]
		if l[4] != '':
			dG.append(float(l[4]))
			seq.append(l[1])

def isref_analysis(opts):
	# reset options list. Fasta will be reset in loop.
	opts['smooth_window'] = 	True
	opts['standard_plot'] = 	False
	opts['auto'] = 			False
	opts['fit_to_sin'] = 		False
	opts['display_sin'] =		False
	opts['r_cutoff'] = 		0.2
	opts['make_palign'] =		False
	opts['print_palign'] = 		False
	opts['use_weblogo'] = 		False
	opts['add_peaks'] = 		False
	opts['print_fasta'] = 		False
	opts['plot_phi_shift'] = 	False
	opts['auto_filter'] = 		True
	opts['ac_cutoff'] = 		-0.2
	opts['plot_auto_filter'] = 	False
	for fil in ['has_motifs.fasta','partial_motifs.fasta','no_motifs.fasta']:
		main.MPEX_tools.reset()
		print "working in", fil, "now!"
		try:
			open(fil,'r')
		except IOError:
			print "Run 'is_reflectin' calculation first. Files do not exist"; exit()
		opts['fasta'] = fil
		with open(fil.split('.')[0]+'_Tic.fasta','w') as f, open(fil.split('.')[0]+'_nTic.fasta','w') as fn:
			for name,seq,dG,i in parse_fasta(opts):
				calc = main_calcs(opts,name,seq,dG,i)
				if calc.period != False:
					f.write(">"+calc.name+"\n"+calc.seq+"\n")
				else:
					fn.write(">"+calc.name+"\n"+calc.seq+"\n")
			if opts['auto_filter']:
				print float(main.MPEX_tools.af)/main.MPEX_tools.nchecked, "% are periodic"
			if display(opts):
				plt.show()
	return None

def main_calcs(opts,name,seq,dG,i):
	calc = main.MPEX_tools(opts,name,seq,dG,i)
	calc.run_calcs()
	return calc

def run(opts):
	if 'Compute' in opts['calcs']:
		for name,seq,dG,i in parse_fasta(opts):
			main_calcs(opts,name,seq,dG,i)
	elif 'Results' in opts['calcs']:
		for name,seq,dG,i in parse_MPEX(opts):
			main_calcs(opts,name,seq,dG,i)
	if display(opts):
		plt.show()
	if opts['isref_analysis']:
		print float(main.MPEX_tools.nref)/main.MPEX_tools.nchecked, "% are reflectin"
		isref_analysis(opts)
	elif opts['auto_filter']:
		print float(main.MPEX_tools.af)/main.MPEX_tools.nchecked, "% are periodic"
	return None
