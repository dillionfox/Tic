from __future__ import division
from __future__ import print_function
import re, os
from tic_utils import calcs
from tic_utils import table
from tic_utils import helpstr
from tic_utils import parser
from tic_utils import ref_tests
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
				print("Bio package not found. \npip install Biopython")
				exit()
		if self.sin_filter or self.display_sin or self.plot_phi_shift:
			try:
				import scipy
			except:
				print("scipy is needed for curve fitting. \npip install scipy")
				print("will continue without sin_filter feature")
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
				print("matplotlib required for plotting. will continue without.")
				self.display = False
				self.display_sin = False
		return None

	def __doc__(self):
		return helpstr.str()

	def __repr__(self):
		return '__repr__ for Compute class in Tic module'

	def __str__(self):
		return '__str__ for Compute class in Tic module'


	def display(self):
		if self.sin_filter or self.display_sin or self.plot_phi_shift or self.standard_plot or self.plot_auto_filter or self.plot_auto_filter_rejected:
			return True
		return False

	def ref_analysis(self):
		# reset options list. Fasta will be reset in loop.
		self.smooth_window = True
		self.auto_filter = True
		self.ac_cutoff = -0.2
		for fil in ['has_motifs.fasta','partial_motifs.fasta','no_motifs.fasta']:
			calcs.MPEX_tools.reset()
			print("working in", fil, "now!")
			try:
				open(fil,'r')
			except IOError:
				print("Run 'is_reflectin' calculation first. Files do not exist"); exit()
			self.fasta = fil
			with open(fil.split('.')[0]+'_Tic.fasta','w') as f, open(fil.split('.')[0]+'_nTic.fasta','w') as fn:
				for name,seq,dG,i in parser.parse_fasta(SeqIO,self.fasta,self.dG_scale):
					calc = self.main_calcs(name,seq,dG,i)
					if calc.period != False:
						f.write(">"+calc.name+"\n"+calc.seq+"\n")
					else:
						fn.write(">"+calc.name+"\n"+calc.seq+"\n")
				if self.auto_filter:
					try:
						print(100*float(calcs.MPEX_tools.af)/calcs.MPEX_tools.nchecked, "% are periodic - auto_filter")
					except ZeroDivisionError:
						print("can't divide by zero")
				if self.sin_filter:
					try:
						print(100*float(calcs.MPEX_tools.sf)/calcs.MPEX_tools.nchecked, "% are periodic - sin_filter")
					except ZeroDivisionError:
						print("can't divide by zero")
				if self.display():
					plt.show()
		return None

	def main_calcs(self,name,seq,dG,i):
		calc = calcs.MPEX_tools(self.__dict__,name,seq,dG,i)
		calc.run_calcs()
		return calc
	
	def run(self):
		if 'Compute' in self.calcs:
			for name,seq,dG,i in parser.parse_fasta(SeqIO,self.fasta,self.dG_scale):
				self.main_calcs(name,seq,dG,i)
		elif 'Results' in self.calcs:
			for name,seq,dG,i in parser.parse_MPEX(self.MPEX_fil):
				self.main_calcs(name,seq,dG,i)
		if self.display():
			if 'DISPLAY' in os.environ:
				plt.show()
			else:
				print("Cannot open display. Saving plot as 'result.png' in working directory")
				plt.savefig('result.png')
		if self.isref_analysis or self.auto_filter or self.sin_filter:
			with open("output.txt", "a") as f:
				if self.isref_analysis:
					f.write(str(self.fasta)+":\t\t"+str(self.dG_scale)+":\t\t"+str(100*float(calcs.MPEX_tools.nref)/calcs.MPEX_tools.nchecked)+"%\n")
					self.ref_analysis()
				elif self.auto_filter:
					f.write(str(self.fasta)+":\t\t"+str(self.dG_scale)+":\t\t"+str(100*float(calcs.MPEX_tools.af)/calcs.MPEX_tools.nchecked)+"%\n")
				elif self.sin_filter:
					f.write(str(self.fasta)+":\t\t"+str(self.dG_scale)+":\t\t"+str(100*float(calcs.MPEX_tools.sf)/calcs.MPEX_tools.nchecked)+"%\n")
		calcs.MPEX_tools.reset()
		return None

if __name__ == "__main__":
	calc = Tic(fasta='GenBank_refs.fasta',standard_plot=True)
	calc.run()
