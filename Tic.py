import re, os
from tic_utils import calcs
from tic_utils import table
from tic_utils import helpstr
from tic_utils.settings import settings
from tic_utils.calcs import MPEX_tools
from tic_utils.calcs import plt

class Error(Exception):
	"""
	Base class to handle custom exceptions

	"""
	pass

class InputError(Error):
	"""
	Custom error messages for mis-handling inputs

	"""
	def __init__(self, expr, message):
		self.expr = expr
		self.message = message

class Compute:

	def __init__(self,**kwargs):
		self.set_kwargs(kwargs)
		self.check_input()
		self.import_libs()

	def set_kwargs(self,kwargs):
		self.__dict__.update((key, value) for key,value in zip(settings.opts.keys(),settings.opts.values()))
		self.__dict__.update((key, value) for key, value in kwargs.items() if key in settings.opts.keys())

	def check_input(self):
		if 'Compute' in self.calcs:
			if self.fasta == '':
				raise InputError('fasta','Default mode is "Compute". Must define a .fasta or specify calcs = ["Results"]')
			if self.dG_scale == 'charge' and self.auto_filter:
				raise InputError('auto_filter',"The 'auto_filter' function does not play well with the 'charge' scale. Consider using 'sin_filter' instead")
		if 'Results' in self.calcs:
			if self.MPEX_fil == '':
				raise InputError('MPEX file','Must define an MPEX Results file to work with "Results" feature')
		return None

	def import_libs(self):
		if 'Compute' in self.calcs:
			try:
				global SeqIO
				from Bio import SeqIO
			except ImportError as error:
				print "Bio package not found. \npip install Biopython"
				exit()
		if self.sin_filter or self.display_sin or self.plot_phi_shift:
			try:
				import scipy
			except:
				print "scipy is needed for curve fitting. \npip install scipy"
				print "will continue without sin_filter feature"
				self.sin_filter = False
				self.display_sin = False
				self.plot_phi_shift = False
			if (self.display_sin and not self.sin_filter) or (self.plot_phi_shift and not self.sin_filter):
				self.sin_filter = True
			return None
		if self.display or self.display_sin or self.plot_phi_shift:
			try:
				import matplotlib
			except:
				print "matplotlib required for plotting. will continue without."
				self.display = False
				self.display_sin = False
		return None

	def __doc__(self):
		return helpstr.str()

	def __repr__(self):
		return '__repr__ for Compute class in Tic module'

	def __str__(self):
		return '__str__ for Compute class in Tic module'

	def read_fasta(self):
		for s in SeqIO.parse(self.fasta,'fasta'):
			yield s

	@staticmethod
	def check_seq(seq):
		for aa in seq:
			if aa not in settings.acceptable_aa: 
				return True
		return False

	def display(self):
		if self.standard_plot:
			return True
		if self.sin_filter or self.display_sin or self.plot_phi_shift:
			return True
		if self.auto_filter:
			return True
		return False
	
	def parse_fasta(self):
		window_size = 19 ; dw = window_size//2 ; 
		for counter,seq in enumerate(self.read_fasta()):
			trim_seq = re.sub('-','', str(seq.seq))
			if self.check_seq(trim_seq): continue
			RFE = [table.table(s,self.dG_scale) for s in trim_seq]
			dG = []
			for ind in range(dw,len(trim_seq)-dw+1):
				s = RFE[ind-dw:ind+dw+1]
				if len(s) != window_size: continue
				dG.append(sum(s))
			yield seq.id,trim_seq,dG,counter
	
	def parse_MPEX(self):
		nseq = 60
		counter = 0
		dG = []
		rseq = []
		window_size = 19
		dw = window_size//2
		phi_list = []
		for it,line in enumerate(open(self.MPEX_fil)):
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
	
	def isref_analysis(self):
		# reset options list. Fasta will be reset in loop.
		self.smooth_window = True
		self.auto_filter = True
		self.ac_cutoff = -0.2
		for fil in ['has_motifs.fasta','partial_motifs.fasta','no_motifs.fasta']:
			calcs.MPEX_tools.reset()
			print "working in", fil, "now!"
			try:
				open(fil,'r')
			except IOError:
				print "Run 'is_reflectin' calculation first. Files do not exist"; exit()
			self.fasta = fil
			with open(fil.split('.')[0]+'_Tic.fasta','w') as f, open(fil.split('.')[0]+'_nTic.fasta','w') as fn:
				for name,seq,dG,i in self.parse_fasta():
					calc = self.main_calcs(name,seq,dG,i)
					if calc.period != False:
						f.write(">"+calc.name+"\n"+calc.seq+"\n")
					else:
						fn.write(">"+calc.name+"\n"+calc.seq+"\n")
				if self.auto_filter:
					print 100*float(calcs.MPEX_tools.af)/calcs.MPEX_tools.nchecked, "% are periodic - auto_filter"
				if self.sin_filter:
					print 100*float(calcs.MPEX_tools.sf)/calcs.MPEX_tools.nchecked, "% are periodic - sin_filter"
				if self.display():
					plt.show()
		return None

	def main_calcs(self,name,seq,dG,i):
		calc = calcs.MPEX_tools(self.__dict__,name,seq,dG,i)
		calc.run_calcs()
		return calc
	
	def run(self):
		if 'Compute' in self.calcs:
			for name,seq,dG,i in self.parse_fasta():
				self.main_calcs(name,seq,dG,i)
		elif 'Results' in self.calcs:
			for name,seq,dG,i in self.parse_MPEX():
				main_calcs(name,seq,dG,i)
		if self.display():
			plt.show()
		if self.isref_analysis:
			print 100*float(calcs.MPEX_tools.nref)/calcs.MPEX_tools.nchecked, "% are reflectin"
			self.isref_analysis()
		elif self.auto_filter:
			print 100*float(calcs.MPEX_tools.af)/calcs.MPEX_tools.nchecked, "% are periodic - auto_filter"
		elif self.sin_filter:
			print 100*float(calcs.MPEX_tools.sf)/calcs.MPEX_tools.nchecked, "% are periodic - sin_filter"
		return None

#if __name__ == "__main__":
#	calc = Tic(fasta='path/to/filename.fasta',standard_plot=True,display=True)
#	calc.run()
