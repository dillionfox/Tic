from utils import Compute as C
from utils import ResultsSummary as RS

calcs = ["Compute", "~Results","~ResultsSummary"] 						# "Results", "ResultsSummary"
# MPEX Results Summary, Results Summary
#MPEX_fil = "MPEX_Results/All_tree_selections-4newgblocks_MPExResultsSummary.txt" ; MPEX_opts = ["~inflate", "print table"]
MPEX_fil = "MPEX_Results/All_tree_selections-4newgblocks_MPExResults.txt" ; MPEX_opts = []
fasta = "/home/dillion/uniprot_sprot.fasta"
weblogo_txt = "ABCDEFG"
ex_seq = "ABCDEFG"

opts = {'calcs':		calcs,\
	'fasta':		fasta,\
	'MPEX_fil':		MPEX_fil,\
	'MPEX_opts':		MPEX_opts,\
	'dG_scale':		0,\
	'smooth_window':	False,\
	'standard_plot':	False,\
	'auto':			False,\
	'is_reflectin':		False,\
	'isref_analysis':	False,\
	'fit_to_sin':		False,\
	'display_sin':		False,\
	'r_cutoff':		0.2,\
	'make_palign':		False,\
	'print_palign':		False,\
	'use_weblogo':		False,\
	'weblogo_txt':		weblogo_txt,\
	'ex_seq':		ex_seq,\
	'add_peaks':		False,\
	'print_fasta':		False,\
	'plot_phi_shift':	False,\
	'auto_filter':		True,\
	'ac_cutoff':		-0.3,\
	'plot_auto_filter':	False}

if __name__ == "__main__":
	if "Compute" in calcs or "Results" in calcs:
		C.run(opts)
	if "ResultsSummary" in calcs:
		seqs = RS.read_ResultsSummary(opts['MPEX_fil'])
		if len(seqs) > 0:
			for s in seqs:
				s.stats()
				if "print table" in opts['MPEX_opts']: s.print_table()
		if "inflate" in opts['MPEX_opts']:
			new_seqs = RS.inflate_seqs(seqs)
			for n in new_seqs:
				print n.name,n.sequence,n.dG,n.rivet,n.RandP,n.aromatics
