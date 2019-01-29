class seq:
	_counter = 0
	def __init__(self,name,sequence=None,dG=None,rivet=None,RandP=None,aromatics=None):
		def input_handler(i):
			if i == None:
				return []
			return i
		# count instances of object
		seq._counter += 1
		self.id = seq._counter
		# general attributes
		self.name = name
		self.sequence = input_handler(sequence)
		self.dG = input_handler(dG)
		self.rivet = input_handler(rivet)
		self.RandP = input_handler(RandP)
		self.aromatics = input_handler(aromatics)
		self.len = []

	def __len__(self):
		return len(self.sequence)

	def __str__(self):
		return "Object containing all computed attributes of subsequence"

	def compute_r(self):
		for s in self.sequence:
			N = s.count("N")
			Q = s.count("Q")
			S = s.count("S")
			total = N+Q+S
			self.rivet.append(total)
		return None

	def compute_len(self):
		for s in self.sequence:
			self.len.append(len(s))
		return None

	def compute_a(self):
		for s in self.sequence:
			Y = s.count("Y")
			F = s.count("F")
			W = s.count("W")
			total = Y+F+W
			self.aromatics.append(total)
		return None

	def compute_RP(self):
		for s in self.sequence:
			R = s.count("R")
			P = s.count("P")
			total = R+P
			self.RandP.append(total)
		return None

	def stats(self):
		for _ in range(len(self.sequence)):
			self.compute_r()
			self.compute_len()
			self.compute_a()
			self.compute_RP()
		return None

	@staticmethod
	def table_header():
		print "name\tdG\trivett\tRandP\taromatics"

	def print_table(self):
		if self.id == 1:
			self.table_header()
		for s in range(len(self.sequence)):
			print self.name,"\t", self.dG[s],"\t", self.rivet[s],"\t", self.RandP[s],"\t", self.aromatics[s]
		return None

