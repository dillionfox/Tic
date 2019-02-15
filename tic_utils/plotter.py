import numpy as np
import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
	mpl.use('Agg')
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
import matplotlib.pyplot as plt
from tic_utils import autotools

def clearplt():
	plt.clf()
	return None

def savefig(name):
	plt.savefig(name)
	return None

def showplt():
	plt.show()
	return None

def display(fasta,scale='corr',name="default"):
	plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9, top=0.9)
	if name == "default":
		name = fasta.split('.')[0]+'_'+scale+'_plot.png'
	if 'DISPLAY' in os.environ: 
		showplt()
	else:
		print("Cannot open display. Saving plot as 'result.png' in working directory")
		savefig(name)
	return None

def plot_shifted(dG,shift):
	x = np.linspace(0,len(dG),num=len(dG))-shift
	plt.ylabel("dG (smoothed)")
	plt.xlabel("Residue Number")
	plt.plot(x,dG)
	return None

def plot(x,auto,label="dG"):
	#---Standard plotting function
	plt.xlabel("Residue Number")
	if not auto:
		if label == "dG":
			plt.ylabel("dG (smoothed)")
		else:
			plt.ylabel(label)
		plt.plot(x)
		return None
	plt.ylabel("Autocorrelation")
	plt.plot(autotools.autocorr(x))
	return None

def sin_plot(x,y,yaj):
	plt.plot(x,y,'--',x,yaj, 'r-')
	return None
