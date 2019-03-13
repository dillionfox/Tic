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

def display(fasta,name="default",save=False):
	plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9, top=0.9)
	if name == "default":
		name = fasta.split('.')[0]+'_'+'_plot.png'
	if 'DISPLAY' in os.environ and save == False: 
		showplt()
	else:
		print("Cannot open display. Saving plot as 'result.png' in working directory")
		savefig(name+'.png')
	return None

def scatter(x,y,c,label,mode=1):
	if mode == 1:
		plt.scatter(x,y,c=c,label=label)
	elif mode == 2:
		plt.scatter(x,y,facecolors='none',s=100,linewidths=0.3,edgecolors='k')
	elif mode == 3:
		plt.scatter(x,y,facecolors='none',s=150,linewidths=0.3,edgecolors='k',marker='s')
	plt.xlabel('Mean')
	plt.ylabel('Variance')
	plt.legend()
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
			#decorate_plot()
			plt.ylabel("dG (smoothed)")
		else:
			plt.ylabel(label)
		plt.plot(x)
		return None
	plt.ylabel("Autocorrelation")
	plt.plot(autotools.corr(x,x))
	return None

def hist(x):
	plt.hist(x,edgecolor='k')

def sin_plot(x,y,yaj):
	plt.plot(x,y,'--',x,yaj, 'r-')
	return None

def decorate_plot():
	import re
	#---No longer being used
	seq = 'MSSNNMWGNMSNNNMWGNMSNNNMWGNNMSGNMFNNNMWGNMNRGRYRGMMEPMSRMTMDFQGRYMDSCGRMVDPRFNDYYGRWNDYDRYYGRSMFNYGWMMNGDRYNRNFRSMDFPERYMDMSGYQMDMCGRWMDPSGRQCNPFNQWSYNRHGCYPGYSYGRNMCYPERWMDMSNYSMDMQGRYMDRWGRQCNPFSQYMNYYGRYWNYPGYNNYYNRNMSYPERHFDMSSWQMDMQGRWMDMQGRYNSPYWSNWYGRNMYNPCQNNQWYGRGDYPGMDCSNWQMDMQGRGMDMQGRGMDMQGRGMDMQGGYMNSWMGDSCYNNW'

	fig,ax = plt.subplots()

	PER = [m.start()-10 for m in re.finditer('PER', seq)]
	mlen = 24

	import matplotlib.patches as patches
	# internal motifs - motif 1
	rect = patches.Rectangle((PER[0],-8),mlen,16,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
	ax.add_patch(rect)
	# linker 1
	rect = patches.Rectangle((PER[0]+mlen,-8),PER[1]-(PER[0]+mlen),16,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
	ax.add_patch(rect)
	# motif 2
	rect = patches.Rectangle((PER[1],-8),mlen,16,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
	ax.add_patch(rect)
	# linker 2
	rect = patches.Rectangle((PER[1]+mlen,-8),PER[2]-(PER[1]+mlen),16,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
	ax.add_patch(rect)
	# motif 3
	rect = patches.Rectangle((PER[2],-8),mlen,16,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
	ax.add_patch(rect)
	# linker 3
	rect = patches.Rectangle((PER[2]+mlen,-8),25,16,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
	ax.add_patch(rect)

	return None
