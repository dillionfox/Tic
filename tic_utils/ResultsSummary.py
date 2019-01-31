import seq

def read_ResultsSummary(fil):
	seqs = []
	for line in open(fil):
		try:
			l = line.split()
			first = l[0]
		except:
			continue
		# only consider lines with things in them
		if len(l) == 0:
			continue
		# extract name for each sequence
		if l[0] == "Results":
			name = l[4]
			seqs.append(seq.seq(name))
		if len(l)<2:
			continue
		if l[1][0] == "(":
			seqs[-1].sequence.append(l[3])
		if "DG" in l and "segment" in l:
			seqs[-1].dG.append(l[2])
	return seqs

def inflate_seqs(seqs):
	new_seqs = []
	for s in seqs:
		for m in range(len(s.sequence)):
			test = seq.seq(s.name,s.sequence[m],s.dG[m],s.rivet[m],s.RandP[m],s.aromatics[m])
			new_seqs.append(test)
	return new_seqs
