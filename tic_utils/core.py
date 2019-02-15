from __future__ import division
from __future__ import print_function
import numpy as np
import os, re
from tic_utils import ref_tests as rt
from tic_utils import autotools
from tic_utils import plotter
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
	mpl.use('Agg')
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
import matplotlib.pyplot as plt

class Tic_tools:

	nchecked = 0	# number of sequences checked
	nref = 0	# number of sequences with motifs
	af = 0		# number of sequences that made it through 'auto_filter'
	sf = 0		# number of sequences that made it through 'sin_filter'

	def __init__(self, opts, name, seq, dG, counter):
		self.opts = opts
		self.name = name
		self.seq = seq
		self.dG = dG
		self.counter = counter
		self.period = None
		self.phi = None
		self.shift = None
		Tic_tools.nchecked += 1

	@property
	def mean(self):
		return np.mean(np.array(self.dG))

	@property
	def var(self):
		return np.var(np.array(self.dG))

	@classmethod
	def reset(cls):
		cls.nchecked = 0
		cls.nref = 0
		cls.af = 0
		cls.sf = 0

	def run_calcs(self):
		if self.opts['is_reflectin']:
			#---Reflectin analysis: Search for motifs.
			Tic_tools.nref+=rt.ref_search(self.opts['fasta'],self.name,self.seq,Tic_tools.nchecked)
		if self.opts['smooth_window']:
			#---Smooth window hydropathies with scaled window smoothing
			self.smooth()
		if self.opts['auto_filter']:
			#---Remove all sequences that don't have at least a 100 residue periodic section
			self.period = autotools.auto_filter(self.dG,self.opts)
			if self.period: Tic_tools.af+=1
			if self.opts['plot_auto_filter'] and self.period != False:
				plotter.plot(self.dG,self.opts['auto'])
			elif self.opts['plot_auto_filter_rejected'] and self.period == False:
				plotter.plot(self.dG,self.opts['auto'])
		if self.opts['sin_filter'] or self.opts['display_sin'] or self.opts['plot_phi_shift']:
			#---Attempt to fit dG to a sine function. Return phase shift
			try:
				self.sin_fit()
			except RuntimeError:
				self.phi = None
			if self.phi != False and self.phi != None:
				if self.opts['print_fasta']: self.make_readable_fasta()
				if self.opts['plot_phi_shift']: plotter.plot_shifted(self.dG, self.phi)
				Tic_tools.sf+=1
		if self.opts['make_palign']:
			#---Make pseudo-alignment of sequences by lining them up by location of 'EPM' residues
			self.pseudo_align()
			if self.shift != False and self.shift != None:
				plotter.plot_shifted(self.dG,-1*self.shift)
		if self.opts['standard_plot']:
			#---plot dG or autocorrelation vs. residue number. Depends on utils/settings
			plotter.plot(self.dG,self.opts['auto'])
		return None

	def smooth(self):
		#---Apply additional smoothing to the already smoothed energy curve
		window_size = 5
		dw = window_size//2
		scaled = []
		for i in range(dw,len(self.dG)-dw+1):
			s = self.dG[i-dw:i+dw+1]
			if len(s) == window_size:
				scaled.append(sum(s)/len(s))
		self.dG = scaled
		return None
	
	def make_readable_fasta(self):
		#---Prints the sequences aligned by location of EPM.
		if self.opts['fasta'] != 'aligned.fasta':
			raise Exception('this method was not designed for this file.')
		pre_insertions = ''.join(['-' for _ in range(35-self.phi)])
		post_insertions = ''.join(['-' for _ in range(389-((35-self.phi)+len(self.seq)))])
		complete_seq = pre_insertions+self.seq+post_insertions
		print(">", self.counter)
		print(complete_seq)
		return None
	
	def pseudo_align(self):
		#---aligning by location of EPM. Max is 50
		MAX_EPM = 70
		try:
			EPM_location = self.seq.index('EPM')
		except ValueError:
			return False
		shift_by = MAX_EPM-EPM_location
		insertions = ''.join(['-' for _ in range(shift_by)])
		if self.opts['print_palign']:
			print(">")
			print(insertions+self.seq)
		self.shift = shift_by
		return None
	
	def sin_fit(self):
		#---Try to fit dG curve to sine function
		import scipy as sy
		from scipy.optimize import curve_fit
		from sklearn.metrics import r2_score
		pts = self.dG
		nbins = range(len(pts))
		x = np.asarray(nbins).ravel()
		y = pts
		def func(x,a,b,c,d): 
			return a*np.sin(b*x+c)+d
		p0 = sy.array([6,1/10.2,0,0])
		try:
			coeffs, matcov = curve_fit(func, x, y, p0)
		except TypeError:
			print("TypeError from curve_fit. Skipping. Usually means not enough points")
			return False
		yaj = func(x, coeffs[0], coeffs[1], coeffs[2],coeffs[3])
		if self.opts['display_sin']:
			plotter.sin_plot(x,y,yaj)
			#plt.plot(x,y,'--',x,yaj, 'r-')
		if r2_score(y, yaj) > self.opts['r_cutoff']:
			phi = -1*coeffs[2]/coeffs[1] 
			if coeffs[0] < 0: phi += np.pi/coeffs[1]
			self.phi = int(phi)
			return None
		else: return False
