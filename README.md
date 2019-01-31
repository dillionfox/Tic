# Tic
Code that searches for proteins with (period)ic whole residue "hydropathies"

USAGE:
I don't expect this to be used by anyone else, ever, so I'll keep this brief.
You either have to specify a fasta file (or similar) that can be read by Bio,
or output from MPEx (http://blanco.biomol.uci.edu/mpex/). MPEx is a great tool
but it requires a java interface that I find to be more trouble than it's worth,
and it's not flexible or extendable. My code can do the things I would otherwise
use MPEx for plus everything I want to do with the output, and it's easy to add
functionality.

*There are lots of dependencies but all are pretty common:*\n
	numpy\n
	matplotlib\n
	os\n
	re\n
	scipy (can be avoided)\n

*There are several ways to analyze output:*\n
	standard_plot	Plot the output (dG vs. residue positions), or the\n
			autocorrelation of dG\n
	auto		Compute autocorrelation of data\n
	smooth_window	Further smooth dG (which is already smoothed)\n
	auto_filter	Throw away sequences that don't look periodic (according\n
			to autocorrelation function)\n
	fit_to_sin	Use SciPy to try to fit a sine wave to dG\n

*All default options are stored in tic_utils/settings.py. To use function, change value to True\n
	***Option		Default Value 	Description**\n
	'calcs':		['Compute']	Internal flag. Can also include\n
						'Results' if analyzing output from MPEx\n
	'fasta':		None		Name of fasta file. Not required if\n
						analyzing MPEx output.\n
	'MPEX_fil':		None		Name of MPEx file. Not required if\n
						analyzing fasta file.\n
	'MPEX_opts':		None		Currently not being used.\n
	'dG_scale':		'interface'	Hydropathy scales: interface, octanol,\n
						octanol-interface, kyte-doolittle, charge\n
	'smooth_window':	False		Smooth the already smoothed curves\n
	'standard_plot':	False		Plot either dG or autocorrelation(dG)\n
	'auto':			False		Choose either dG or auto(dG)\n
	'fit_to_sin':		False		Fit dG curve to sine wave\n
	'display_sin':		False		Display curve fits (unaligned)\n
	'plot_phi_shift':	False		Align fitted curves by phase shift (and plot)\n
	'r_cutoff':		0.2		r^2 cutoff to accept fit\n
	'use_weblogo':		False		Print a consensus sequence as x-axis of plot\n
	'weblogo_txt':		''		Text\n
	'add_peaks':		False		Use SciPy function to find peaks in autocorr plot\n
	'auto_filter':		False		Filter out sequences with aperiodic dG curves\n
	'ac_cutoff':		-0.3		Cutoff for auto_filter. Lower numbers = more periodic\n
	'plot_auto_filter':	False		Show periodic sequences only\n
	'is_reflectin':		False		Special case\n
	'isref_analysis':	False		Special case\n
	'make_palign':		False		Special case\n
	'print_palign':		False		Special case\n
	'ex_seq':		''		Special case\n
	'print_fasta':		False		Special case. Print fasta containing sequences\n
