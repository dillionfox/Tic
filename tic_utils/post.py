from __future__ import print_function
from tic_utils.core import Tic_tools
from tic_utils import autotools
from tic_utils import plotter
import numpy as np

class PostProc:
	def __init__(self,data):
		self.data = np.array(data)
		self.corr = None
		plotter.clearplt()

	def print_data(self):
		print("placeholder function. This code is under development")
		# example usage of correlation function between two data sets
		#c = autotools.corr(self.data[0][0],self.data[0][1])
		#print(c)

		# example of how to access data
		#for fil in self.data:
		#	for dG in fil:
		#		print(dG.var)

	def run_corr(self):
		for i in range(len(self.data[0])):
			self.corr = autotools.corr(self.data[0][i].dG, self.data[1][i].dG)
			plotter.plot(self.corr,False,'Correlation')
		return None

	def run(self):
		#self.print_data()
		self.run_corr()
		return None

	def plot(self):
		plotter.plot(self.corr,False)
