# Tic
Code that searches for proteins with (period)ic whole residue "hydropathies"

USAGE:
I don't expect this to be used by anyone else, ever, so I'll keep this brief.

	import Tic

	calc = Tic(fasta="/path/to/file.fasta",standard_plot=True,...)

	calc.run()

You either have to specify a fasta file (or similar) that can be read by Bio,
or output from MPEx (http://blanco.biomol.uci.edu/mpex/). MPEx is a great tool
but it requires a java interface that I find to be more trouble than it's worth,
and it's not flexible or extendable. My code can do the things I would otherwise
use MPEx for plus everything I want to do with the output, and it's easy to add
functionality.

There are lots of dependencies but all are pretty common:

	Python 2.7

	numpy

	matplotlib

	os

	re

	scipy (can be avoided)

There are several ways to analyze output:

	standard_plot	Plot the output (dG vs. residue positions), or the
			autocorrelation of dG

	auto		Compute autocorrelation of data

	smooth_window	Further smooth dG (which is already smoothed)

	auto_filter	Throw away sequences that don't look periodic (according
			to autocorrelation function)

	fit_to_sin	Use SciPy to try to fit a sine wave to dG

*All default options are stored in tic_utils/settings.py. To use function, change value to True

	Option			Default Value 	Description

	'calcs':		['Compute']	Internal flag. Can also include
						'Results' if analyzing output from MPEx

	'fasta':		None		Name of fasta file. Not required if
						analyzing MPEx output.

	'MPEX_fil':		None		Name of MPEx file. Not required if
						analyzing fasta file.

	'MPEX_opts':		None		Currently not being used.

	'dG_scale':		'interface'	Hydropathy scales: interface, octanol,
						octanol-interface, kyte-doolittle, charge

	'smooth_window':	False		Smooth the already smoothed curves

	'standard_plot':	False		Plot either dG or autocorrelation(dG)

	'auto':			False		Choose either dG or auto(dG)

	'fit_to_sin':		False		Fit dG curve to sine wave

	'display_sin':		False		Display curve fits (unaligned)

	'plot_phi_shift':	False		Align fitted curves by phase shift (and plot)

	'r_cutoff':		0.2		r^2 cutoff to accept fit

	'use_weblogo':		False		Print a consensus sequence as x-axis of plot

	'weblogo_txt':		''		Text

	'add_peaks':		False		Use SciPy function to find peaks in autocorr plot

	'auto_filter':		False		Filter out sequences with aperiodic dG curves

	'ac_cutoff':		-0.3		Cutoff for auto_filter. Lower numbers = more periodic

	'plot_auto_filter':	False		Show periodic sequences only

	'is_reflectin':		False		Special case

	'isref_analysis':	False		Special case

	'make_palign':		False		Special case

	'print_palign':		False		Special case

	'ex_seq':		''		Special case

	'print_fasta':		False		Special case. Print fasta containing sequences
