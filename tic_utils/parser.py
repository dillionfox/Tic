from tic_utils.settings import settings
import table

def read_fasta(SeqIO,fasta):
	for s in SeqIO.parse(fasta,'fasta'):
		yield s

def parse_fasta(SeqIO,fasta,dG_scale):
	import re
	window_size = 19 ; dw = window_size//2 ; 
	for counter,seq in enumerate(read_fasta(SeqIO,fasta)):
		trim_seq = re.sub('-','', str(seq.seq))
		if check_seq(trim_seq): continue
		if settings.opts['skip_100'] and len(trim_seq) > 100:
			trim_seq = trim_seq[100:]
		RFE = [table.table(s,dG_scale) for s in trim_seq]
		dG = []
		for ind in range(dw,len(trim_seq)-dw+1):
			s = RFE[ind-dw:ind+dw+1]
			if len(s) != window_size: continue
			dG.append(sum(s))
		yield seq.id,trim_seq,dG,counter

def parse_MPEX(MPEX_fil):
	nseq = 60
	counter = 0
	dG = []
	rseq = []
	window_size = 19
	dw = window_size//2
	phi_list = []
	for it,line in enumerate(open(MPEX_fil)):
		if counter >= nseq: break
		try:
			l = line.split('\t') ; first = l[0]
		except:
			continue
		if len(l) == 0: continue
		if l[0] == '"Position"':
			if len(dG) > 0:
				counter += 1
				fullseq = ''.join(seq).upper() 
				yield it,fullseq,dG,it
			dG = [] ; seq = [] ; rseq = [] ; continue
		resname = l[1]
		if l[4] != '':
			dG.append(float(l[4]))
			seq.append(l[1])

def check_seq(seq):
	for aa in seq:
		if aa not in settings.acceptable_aa: 
			return True
	return False
