from __future__ import print_function
from tic_utils.core import Tic_tools
from tic_utils import autotools
from tic_utils import plotter
import numpy as np

class PostProc:
	def __init__(self,data,task,fasta,opts):
		self.data = np.array(data)
		self.corr = None
		self.task = task
		self.fasta = fasta
		self.opts = opts
		plotter.clearplt()
		#self.check_input()

	def check_input(self):
		if self.data.shape[0] != 2 and "correlate" in self.task:
			print("You must have two dG scales to run correlation")
			exit()
		return None

	def run_corr(self):
		for i in range(len(self.data[0])):
			self.corr = autotools.corr(self.data[0][i].dG, self.data[1][i].dG)
			plotter.plot(self.corr,False,'Correlation')
		return None

	def corr_2(self):
		if self.data.shape[0] == 1 and self.data.shape[1] == 2:
			self.corr = autotools.corr(self.data[0][0].dG, self.data[0][1].dG)
			plotter.plot(self.corr,False,'Correlation')

	def mean_var_scatter(self):
		# also circle charge ones that are periodic in interface and vice versa
		m = [dat.mean for dat in self.data[0]]
		v = [dat.var  for dat in self.data[0]]
		plotter.scatter(m,v,'b',"Interfelicity")

		m = [dat.mean for dat in self.data[0] if dat.period == True]
		v = [dat.var  for dat in self.data[0] if dat.period == True]
		plotter.scatter(m,v,'b',"Interfelicity",2)

		m = [dat.mean for i,dat in enumerate(self.data[0]) if self.data[1][i].period == True]
		v = [dat.var  for i,dat in enumerate(self.data[0]) if self.data[1][i].period == True]
		plotter.scatter(m,v,'b',"Interfelicity",3)

		m = [dat.mean for dat in self.data[1]]
		v = [dat.var  for dat in self.data[1]]
		plotter.scatter(m,v,'r',"Charge")

		m = [dat.mean for dat in self.data[1] if dat.period == True]
		v = [dat.var  for dat in self.data[1] if dat.period == True]
		plotter.scatter(m,v,'r',"Charge",2)

		m = [dat.mean for i,dat in enumerate(self.data[1]) if self.data[0][i].period == True]
		v = [dat.var  for i,dat in enumerate(self.data[1]) if self.data[0][i].period == True]
		plotter.scatter(m,v,'b',"Charge",3)

		plotter.display(self.fasta,'Mean_vs_Var',self.opts['save'])
		return None

	def write_all(self):
		[f1,f2,f3,f4] = [self.fasta.split('.')[0]+i for i in ['_both.fasta','_interface.fasta','_charge.fasta','_neither.fasta']]
		with open(f1,'w') as b1, open(f2,'w') as b2, open(f3,'w') as b3, open(f4,'w') as b4 :
			for i in range(self.data.shape[1]):
				bool1 = self.data[0][i].period ; bool2 = self.data[1][i].period ; out = ">"+self.data[0][i].name+"\n"+self.data[0][i].seq
				if bool1 and bool2: b1.write(out)
				elif bool1 and not bool2: b2.write(out)
				elif not bool1 and bool2: b3.write(out)
				elif not bool1 and not bool2: b4.write(out)
		return None

	def stats(self):
		plotter.clearplt()
		self.mean_var_scatter()
		self.write_all()
		return None

	def run(self):
		if "correlate" in self.task:
			#self.run_corr()
			self.stats()
		if "corr_2" in self.task:
			self.corr_2()
		return None

	def plot(self):
		plotter.plot(self.corr,False)
