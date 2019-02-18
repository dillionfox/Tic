from __future__ import print_function
from tic_utils.core import Tic_tools
from tic_utils import autotools
from tic_utils import plotter
import numpy as np

class PostProc:
	def __init__(self,data,task):
		self.data = np.array(data)
		self.corr = None
		self.task = task
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

	def stats(self):
		for i in range(len(self.data)):
			print(i, self.data[i].shape)

	def run(self):
		if "correlate" in self.task:
			self.run_corr()
			self.stats()
		if "corr_2" in self.task:
			self.corr_2()
		return None

	def plot(self):
		plotter.plot(self.corr,False)
