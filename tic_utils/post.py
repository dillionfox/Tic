from __future__ import print_function
from tic_utils.core import Tic_tools
from tic_utils import autotools
import numpy as np

class PostProc:
	def __init__(self,data):
		self.data = data

	def print_data(self):
		print("placeholder function. Under development")
		# example usage of correlation function between two data sets
		#c = autotools.corr(self.data[0][0],self.data[0][1])
		#print(c)

		# example of how to access data
		for fil in self.data:
			for dG in fil:
				print(dG.var)

	def run(self):
		self.print_data()
