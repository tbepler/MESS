import math
from collections import Counter


class PositionWeightMatrix:
	def __init__(self,scores,length):
		self._struct = scores
		self._n = length
		self._alphabet = scores.keys()
		self._alphabet.sort()
	def __eq__(self,other):
		return pwmEq(self._struct, other._struct, 0)
	def __hash__(self):
		h = 5
		for c in self._alphabet:
			for j in range(len(self._struct[i])):
				h = h*11 + hash(self._struct[c][j])
		return h
	def __getitem__(self, i):
		return self._struct[i]
	def length(self):
		return self._n
	def __len__(self):
		return self._n
	def score(self, s):
		assert len(s) == self.length()
		return sum( [ self._struct[s[i]][i] for i in range(len(s)) ] )
	def scoreAll(self, s):
		assert len(s) >= self.length()
		return [ self.score( s[i:i+self.length()] ) for i in range( len(s) - self.length() + 1 ) ]
			
	def alphabet(self):
		return self._alphabet
	def adjustToBG( self, bg ):
		return adjustToBG( self, bg )
	def __str__(self):
		L = self.length()
		s = "  "
		for i in range(L):
			s += "%-5d " % (i+1)
		s += "\n"
		for c in self._alphabet:
			s += " %s" % c
			for i in range(L):
				s += " %5.3f" % self._struct[c][i]
			s += "\n"
		return s

def pwmEq(pwm1,pwm2,tolerance):
    for i in pwm1.keys():
        if len(pwm1[i]) == len(pwm2[i]):
            for j in range(len(pwm1[i])):
                if abs(exp(pwm1[i][j]) - exp(pwm2[i][j])) > tolerance: return False
        else:
            return False
    return True

def adjustToBG(pwm,bg):
	return PositionWeightMatrix(
		{ c : [ pwm[c][j] - log(bg[c]) for j in range(len(pwm)) ] for c in pwm.alphabet() },
		pwm.length() 
		)

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

def nucleotideFrequency( seq, pseudocount = 0 ):
	counts = Counter( seq.upper().strip() )
	addedCounts = len( counts ) * pseudocount;
	return { x : float( counts[x] + pseudocount ) / float( len( seq ) + addedCounts ) for x in counts.keys() }

def parseProbabilities( f ):
	m = {}
	l = -1
	for line in f:
		tokens = line.strip().split()
		scores = [ log( float( tokens[i] ) ) for i in range( 1, len(tokens ) ) ]
		if l < 0:
			l = len( scores )
		else:
			assert l == len( scores )
		m[ tokens[0][0] ] = scores
	return PositionWeightMatrix( m, l )

def parsePFM( f, pseudocount = 0 ):
	bases = f.readline().strip().split()
	pfm = [ { bases[i]:(float(s)+pseudocount) for (i,s) in enumerate(line.strip().split()) } for line in f ]
	l = len(pfm)
	pfm = { c:[ log( pfm[j][c] ) - log( float(sum([ pfm[j][k] for k in pfm[j].keys() ])) ) for j in range(len(pfm)) ] for c in bases }
	return PositionWeightMatrix(pfm,l)
