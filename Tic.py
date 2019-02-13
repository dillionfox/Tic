from __future__ import division
from __future__ import print_function
import re, os
from Bio import SeqIO
from tic_utils import calcs
from tic_utils import table
from tic_utils import helpstr
from tic_utils import parser
from tic_utils import ref_tests
from tic_utils.settings import settings
from tic_utils.calcs import MPEX_tools
from tic_utils.calcs import plt
from tic_utils.err import InputError

class Compute:

	def __init__(self,**kwargs):
		self.set_kwargs(kwargs)
		self.check_input()

	def set_kwargs(self,kwargs):
		self.__dict__.update((key, value) for key,value in zip(settings.opts.keys(),settings.opts.values()))
		self.__dict__.update((key, value) for key, value in kwargs.items() if key in settings.opts.keys())

	def __doc__(self):
		return helpstr.str()

	def __repr__(self):
		return '__repr__ for Compute class in Tic module'

	def __str__(self):
		return '__str__ for Compute class in Tic module'

	def check_input(self):
		if 'Compute' in self.calcs:
			if self.fasta == '':
				raise InputError('fasta','Default mode is "Compute". Must define a .fasta or specify calcs = ["Results"]')
		if 'Results' in self.calcs:
			if self.MPEX_fil == '':
				raise InputError('MPEX file','Must define an MPEX Results file to work with "Results" feature')
		return None

	def display(self):
		if self.sin_filter or self.display_sin or self.plot_phi_shift or self.standard_plot or self.plot_auto_filter or self.plot_auto_filter_rejected:
			if 'DISPLAY' in os.environ:
				plt.show()
			else:
				print("Cannot open display. Saving plot as 'result.png' in working directory")
				plt.savefig(self.fasta.split('.')[0]+'_'+self.dG_scale+'_plot.png')
		return None

	def print_stats(self):
		def save(x):
			with open("output.txt", "a") as f:
				try:
					f.write(str(self.fasta)+":\t\t"+str(self.dG_scale)+":\t\t"+str(100*x/calcs.MPEX_tools.nchecked)+"%\n")
				except ZeroDivisionError:
					print("can't divide by zero")
		if self.isref_analysis: save(calcs.MPEX_tools.nref)
		if self.auto_filter: save(calcs.MPEX_tools.af)
		if self.sin_filter: save(calcs.MPEX_tools.sf)
		return None

	def ref_analysis(self):
		self.smooth_window = True ; self.auto_filter = True
		for fil in ['has_motifs.fasta','partial_motifs.fasta','no_motifs.fasta']:
			calcs.MPEX_tools.reset()
			print("working on", fil, "now!")
			try:
				open(fil,'r')
			except IOError:
				print("Run 'is_reflectin' calculation first. Files do not exist"); exit()
			self.fasta = fil
			with open(fil.split('.')[0]+'_Tic.fasta','w') as fp, open(fil.split('.')[0]+'_nTic.fasta','w') as fn:
				for name,seq,dG,i in parser.parse_fasta(SeqIO,self.fasta,self.dG_scale):
					calc = self.main_calcs(name,seq,dG,i)
					if calc.period != False:
						fp.write(">"+calc.name+"\n"+calc.seq+"\n")
					else:
						fn.write(">"+calc.name+"\n"+calc.seq+"\n")
				self.print_stats()
				self.display()
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
		self.display()
		if self.isref_analysis: self.ref_analysis()
		self.print_stats()
		calcs.MPEX_tools.reset()
		return None

if __name__ == "__main__":
	calc = Compute(fasta='GenBank_refs.fasta',standard_plot=True)
	calc.run()
