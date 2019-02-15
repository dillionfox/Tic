from __future__ import division
from __future__ import print_function
import re, os
from Bio import SeqIO
import numpy as np
from tic_utils import plotter
from tic_utils import table
from tic_utils import helpstr
from tic_utils import parser
from tic_utils import ref_tests
from tic_utils import core
from tic_utils.core import Tic_tools
from tic_utils.err import InputError
from tic_utils.post import PostProc
from tic_utils.settings import settings

class Tic:
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
		if isinstance(self.fasta,(str,)): self.fasta_list = [self.fasta]
		elif isinstance(self.fasta,(list,)): self.fasta_list = self.fasta
		else: raise InputError('fasta','Must supply a .fasta file')
		for fil in self.fasta_list:
			assert(os.path.isfile(fil))
		if isinstance(self.dG_scale,(str,)): self.dG_scale_list = [self.dG_scale]
		elif isinstance(self.dG_scale,(list,)): self.dG_scale_list = self.dG_scale
		else: raise InputError('dG_scale','See README or settings for possible dG_scale values')
		return None

	def display(self):
		if self.sin_filter or self.display_sin or self.plot_phi_shift or self.standard_plot or self.plot_auto_filter or self.plot_auto_filter_rejected or self.make_palign:
			plotter.display(self.fasta,self.dG_scale)
		self.print_stats()
		return None

	def print_stats(self):
		def save(x):
			with open("output.txt", "a") as f:
				try:
					f.write(str(self.fasta)+":\t\t"+str(self.dG_scale)+":\t\t"+str(100*x/core.Tic_tools.nchecked)+"%\n")
				except ZeroDivisionError:
					print("can't divide by zero")
		if self.isref_analysis: save(core.Tic_tools.nref)
		if self.auto_filter: save(core.Tic_tools.af)
		if self.sin_filter: save(core.Tic_tools.sf)
		return None

	def main_calcs(self,name,seq,dG,i):
		inst = core.Tic_tools(self.__dict__,name,seq,dG,i)
		inst.run_calcs()
		return inst

	def run(self):
		data = []
		for self.fasta in self.fasta_list:
			for self.dG_scale in self.dG_scale_list:
				inst_data = []
				core.Tic_tools.reset() ; plotter.clearplt()
				with open(self.fasta.split('.')[0]+'_Tic.fasta','w') as fp, open(self.fasta.split('.')[0]+'_nTic.fasta','w') as fn:
					for name,seq,dG,i in parser.parse_fasta(SeqIO,self.fasta,self.dG_scale):
						calc = self.main_calcs(name,seq,dG,i)
						if calc.period != False: fp.write(">"+name+"\n"+seq+"\n")
						else: fn.write(">"+name+"\n"+seq+"\n")
						inst_data.append(calc)
				self.display()
				data.append(inst_data)
		post = PostProc(data)
		post.run()
		#post.plot()
		plotter.display(self.fasta,'corr','Correlation')
		return None

if __name__ == "__main__":
	tic = Tic(fasta='GenBank_refs.fasta',standard_plot=True)
	tic.run()
