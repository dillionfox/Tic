from tic_utils.settings import settings
from tic_utils import table

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

def check_seq(seq):
	for aa in seq:
		if aa not in settings.acceptable_aa: 
			return True
	return False
