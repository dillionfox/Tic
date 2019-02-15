import numpy as np
import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
	mpl.use('Agg')
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
import matplotlib.pyplot as plt
from tic_utils import autotools
	
def plot_shifted(dG,shift):
	x = np.linspace(0,len(dG),num=len(dG))-shift
	plt.ylabel("dG (smoothed)")
	plt.xlabel("Residue Number")
	plt.plot(x,dG)
	return None

def plot(dG,auto):
	#---Standard plotting function
	plt.ylabel("dG (smoothed)")
	plt.xlabel("Residue Number")
	if not auto:
		plt.plot(dG)
		return None
	plt.plot(autotools.autocorr(dG))
	return None

def sin_plot(x,y,yaj):
	plt.plot(x,y,'--',x,yaj, 'r-')
	return None
