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
from tic_utils.core import Tic_tools
from tic_utils.core import plt
from tic_utils.err import InputError

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
		if self.fasta == '':
			raise InputError('fasta','Must supply a .fasta file')
		if isinstance(self.fasta,(str,)): 
			self.fasta_list = [self.fasta]
		elif isinstance(self.fasta,(list,)): 
			self.fasta_list = self.fasta
		for fil in self.fasta_list:
			assert(os.path.isfile(fil))
		return None

	def display(self):
		if self.sin_filter or self.display_sin or self.plot_phi_shift or self.standard_plot or self.plot_auto_filter or self.plot_auto_filter_rejected:
			if 'DISPLAY' in os.environ:
				plt.show()
			else:
				print("Cannot open display. Saving plot as 'result.png' in working directory")
				plt.savefig(self.fasta.split('.')[0]+'_'+self.dG_scale+'_plot.png')
		self.print_stats()
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

	def main_calcs(self,name,seq,dG,i):
		inst = core.Tic_tools(self.__dict__,name,seq,dG,i)
		inst.run_calcs()
		return calc

	def run(self):
		for self.fasta in self.fasta_list:
			calcs.MPEX_tools.reset()
			with open(self.fasta.split('.')[0]+'_Tic.fasta','w') as fp, open(self.fasta.split('.')[0]+'_nTic.fasta','w') as fn:
				for name,seq,dG,i in parser.parse_fasta(SeqIO,self.fasta,self.dG_scale):
					calc = self.main_calcs(name,seq,dG,i)
					if calc.period != False:
						fp.write(">"+name+"\n"+seq+"\n")
					else:
						fn.write(">"+name+"\n"+seq+"\n")
			self.display()
		return None

if __name__ == "__main__":
	tic = Tic(fasta='GenBank_refs.fasta',standard_plot=True)
	tic.run()
