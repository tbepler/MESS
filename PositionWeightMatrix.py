import math


class PositionWeightMatrix:
	def __init__(self,l):
		self.struct = l
	def asList(self):
		return self.struct
	def __eq__(self,other):
		return pwmEq(self.struct, other.struct, 0)
	def __hash__(self):
		h = 5
		for i in range(len(alphabet)):
			for j in range(len(self.struct[i])):
				h = h*11 + hash(self.struct[i][j])
		return h
	def __getitem__(self, i):
		return self.struct[i]
	def length(self):
		return len(self.struct[0])
	def score(self, s):
		assert len(s) == self.length()
		return sum( [ self.struct[baseIndex[s[i]]][i] for i in range(len(s)) ] )
	def __str__(self):
		L = len(self.struct[0])
		s = "  "
		for i in range(L):
			s += "%-5d " % (i+1)
		s += "\n"
		for j in range(len(alphabet)):
			s += " %s" % alphabet[j]
			for i in range(L):
				s += " %5.3f" % self.struct[j][i]
			s += "\n"
		return s

def adjustToBG(pwm,bg):
	return PositionWeightMatrix([ [ log(pwm[i][j]/bg[i]) for j in range(len(pwm[i])) ] for i in range(len(alphabet))])

def log(x):
	if x == 0: return float("-inf")
	return math.log(x)

def parseBG( bgF ):
	length = 0
	counts = Counter()
	for line in bgF:
		if line[0] != ">":
			length += len(line)
			counts.update(line.upper().strip())
	return [ float(counts[x])/float(length) for x in alphabet ]

def nucleotideFrequency( seq ):
	counts = Counter( seq.upper().strip() )
	

def parsePFM( f, pseudocount = 0 ):
	bases = f.readline().strip().split()
	pfm = [ { bases[i]:(float(s)+pseudocount) for (i,s) in enumerate(line.strip().split()) } for line in f ]
	pfm = [ [ pfm[j][alphabet[i]]/float(sum([ pfm[j][k] for k in pfm[j].keys() ])) for j in range(len(pfm)) ] for i in range(len(alphabet)) ]
	return PositionWeightMatrix(pfm)
