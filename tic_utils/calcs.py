import numpy as np
import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
	mpl.use('Agg')
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
import matplotlib.pyplot as plt
import re
import ref_tests as rt

class MPEX_tools:

	nchecked = 0	# number of sequences checked
	nref = 0	# number of sequences with motifs
	af = 0		# number of sequences that made it through 'autofilter'

	def __init__(self, opts, name, seq, dG, counter):
		self.opts = opts
		self.name = name
		self.seq = seq
		self.dG = dG
		self.counter = counter
		self.period = None
		self.phi = None
		self.shift = None
		self.has_Nterm = False 
		self.has_Intern = False 
		if self.opts['use_weblogo']:
			self.fig,self.ax = plt.subplots(1)
		MPEX_tools.nchecked += 1

	@classmethod
	def reset(cls):
		cls.nchecked = 0
		cls.nref = 0    
		cls.af = 0	    

	def run_calcs(self):
		if self.opts['is_reflectin']:
			"""
			Reflectin analysis: Search for motifs.
	
			"""
			self.ref_search()
		if self.opts['smooth_window']:
			"""
			Smooth window hydropathies with scaled window smoothing
	
			"""
			self.smooth()
		if self.opts['auto_filter']:
			"""
			Remove all sequences that don't have at least a 100 residue periodic section
	
			"""
			self.auto_filter()
			if self.period == False:
				return False
			else:
				MPEX_tools.af+=1
		if self.opts['fit_to_sin']:
			"""
			Attempt to fit dG to a sine function. Return phase shift
	
			"""
			try:
				self.sin_fit()
			except RuntimeError:
				self.phi = None
			if self.phi != False and self.phi != None:
				if self.opts['print_fasta']:
					self.make_readable_fasta()
				if self.opts['plot_phi_shift']:
					self.plot_shifted()
		if self.opts['make_palign']:
			"""
			Make pseudo-alignment of sequences by lining them up by location of 'EPM' residues
	
			"""
			self.pseudo_align()
			if self.shift != False and self.shift != None:
				self.plot_shifted(-1*self.shift)
		if self.opts['standard_plot']:
			"""
			plot dG or autocorrelation vs. residue number. Depends on utils/settings
	
			"""
			self.plot()
		return None
	
	def ref_search(self):
		"""
		Search sequences for characteristic motifs and collect statistics

		"""
		self.has_Nterm = rt.has_N_term(self.seq)
		self.has_Intern = rt.has_motif(self.seq)
		if MPEX_tools.nchecked == 1:
			open('has_motifs.fasta','w').close()
			open('partial_motifs.fasta','w').close()
			open('no_motifs.fasta','w').close()
		elif self.has_Nterm and self.has_Intern:
			MPEX_tools.nref += 1
			with open('has_motifs.fasta','a') as f:
				f.write(">"+self.name+"\n"+self.seq+"\n")
		elif not self.has_Nterm and not self.has_Intern:
			with open('no_motifs.fasta','a') as f:
				f.write(">"+self.name+"\n"+self.seq+"\n")
		else:
			with open('partial_motifs.fasta','a') as f:
				f.write(">"+self.name+"\n"+self.seq+"\n")
		return None

	def smooth(self):
		"""
		Apply additional smoothing to the already smoothed energy curve

		"""
		window_size = 5
		dw = window_size//2
		scaled = []
		for i in range(dw,len(self.dG)-dw+1):
			s = self.dG[i-dw:i+dw+1]
			if len(s) == window_size:
				scaled.append(sum(s)/len(s))
		self.dG = scaled
		return None
	
	def auto_filter(self):
		"""
		Compute autocorrelation for energy function and determine if result
		is 'periodic enough'

		"""
		exclude = 0
		ac = self.autocorr()
		if any(ac) == False:
			self.period = False
			return False
		if min(ac) > self.opts['ac_cutoff']:
			self.period = False
			#return False
			ac,exclude = self.auto_scan()
			if any(ac) == False: return False
		peaks = self.find_peaks(ac)
		if self.opts['plot_auto_filter']:
			heights = [ac[i] for i in peaks]
			plt.scatter(peaks, heights)
			plt.plot(ac)
		if len(peaks) > 0:
			self.period = peaks[0]
		else:
			self.period = False
	
	def auto_scan(self):
		"""
		If autofilter fails, search sequence for periodic region

		"""
		dx = 20
		for i in range(len(self.dG)//dx-4):
			ac = self.autocorr(self.dG[dx*i:])
			if min(ac) <= self.opts['ac_cutoff']:
				return ac,dx*i
		return [False],False
	
	def make_readable_fasta(self):
		"""
		Prints the sequences aligned by location of EPM.

		"""
		if self.opts['fasta'] != aligned.fasta:
			raise Exception('this method was not designed for this file.')
		pre_insertions = ''.join(['-' for _ in range(35-self.phi)])
		post_insertions = ''.join(['-' for _ in range(389-((35-self.phi)+len(self.seq)))])
		complete_seq = pre_insertions+self.seq+post_insertions
		print ">", self.counter
		print complete_seq
		return None
	
	def pseudo_align(self):
		"""
		aligning by location of EPM. Max is 50
		"""
		MAX_EPM = 70
		try:
			EPM_location = self.seq.index('EPM')
		except ValueError:
			return False
		shift_by = MAX_EPM-EPM_location
		insertions = ''.join(['-' for _ in range(shift_by)])
		if self.opts['print_palign']:
			print ">"
			print insertions+self.seq
		self.shift = shift_by
		return None
	
	def sin_fit(self):
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
		coeffs, matcov = curve_fit(func, x, y, p0)
		yaj = func(x, coeffs[0], coeffs[1], coeffs[2],coeffs[3])
		if self.opts['display_sin']:
			plt.plot(x,y,'--',x,yaj, 'r-')
		if r2_score(y, yaj) > self.opts['r_cutoff']:
			phi = -1*coeffs[2]/coeffs[1] 
			if coeffs[0] < 0: phi += np.pi/coeffs[1]
			self.phi = int(phi)
			return None
		else: return False
	
	def autocorr(self, x=[]):
		if len(x) == 0: x = self.dG
		try:
			result = np.correlate(x,x,mode='full')
		except ValueError:
			return [False]
		return result[result.size/2:]/max(result)
	
	def plot_shifted(self,shift=None):
		if shift == None:
			if self.phi != False:
				shift = self.phi
			else:
				return None
		x = np.linspace(0,len(self.dG),num=len(self.dG))-shift
		plt.ylabel("dG (smoothed)")
		plt.xlabel("Residue Number")
		plt.plot(x,self.dG)
		return None
	
	def find_peaks(self, pts):
		import scipy.signal as sig
		x = np.linspace(0,len(pts))
		if self.opts['auto']:
			ac = np.correlate(pts,pts,mode='full')
			pts = ac[len(ac)/2:]/max(ac)
		peaks = sig.find_peaks(pts,height=0,distance=45)[0]
		return peaks
	
	def plot(self):
		plt.ylabel("dG (smoothed)")
		plt.xlabel("Residue Number")
		if not self.opts['auto']:
			plt.plot(self.dG)
			return None
		plt.plot(self.autocorr(self.dG))
		return None
	
	def decorate_plot(self):
		EPM = opts['weblogo_txt'].index('EPM')-1
		PER = [m.start() for m in re.finditer('PER', opts['weblogo_txt'])]
		mlen = 24
	
		import matplotlib.patches as patches
		# N-terminal motif
		rect = patches.Rectangle((EPM,-15),28,30,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
		self.ax.add_patch(rect)
		# N-linker
		rect = patches.Rectangle((EPM+28,-15),PER[0]-(EPM+28),30,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
		self.ax.add_patch(rect)
		# internal motifs - motif 1
		rect = patches.Rectangle((PER[0],-15),mlen,30,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
		self.ax.add_patch(rect)
		# linker 1
		rect = patches.Rectangle((PER[0]+mlen,-15),PER[1]-(PER[0]+mlen),30,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
		self.ax.add_patch(rect)
		# motif 2
		rect = patches.Rectangle((PER[1],-15),mlen,30,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
		self.ax.add_patch(rect)
		# linker 2
		rect = patches.Rectangle((PER[1]+mlen,-15),PER[2]-(PER[1]+mlen),30,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
		self.ax.add_patch(rect)
		# motif 3
		rect = patches.Rectangle((PER[2],-15),mlen,30,linewidth=1,edgecolor='r',facecolor='r',alpha=0.3)
		self.ax.add_patch(rect)
		# linker 3
		rect = patches.Rectangle((PER[2]+mlen,-15),25,30,linewidth=1,edgecolor='r',facecolor='b',alpha=0.3)
		self.ax.add_patch(rect)
		# motifs I used for simulations
		rect = patches.Rectangle((105,-10.95),34,21.15,linewidth=1,edgecolor='k',facecolor='None',linestyle='--')
		self.ax.add_patch(rect)
		rect = patches.Rectangle((103,-10.95),63,21.15,linewidth=1,edgecolor='g',facecolor='None',linestyle='--')
		self.ax.add_patch(rect)
	
		linker_1 = [PER[0]+mlen, PER[1]]
		linker_2 = [PER[1]+mlen, PER[2]]
		linker_3 = [PER[2]+mlen, PER[2]+mlen+25]
		
		# linkers
		for l in [linker_1, linker_2, linker_3]:
			lseq = opts['ex_seq'][l[0]:l[1]]
			print l, lseq
	
		#plt.xticks(range(len(weblogo_txt)), list(weblogo_txt))
		return None
