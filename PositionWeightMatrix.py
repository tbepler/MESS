import math
import random
import multiprocessing
import itertools
from collections import Counter


class PositionWeightMatrix:

	def __init__(self,scores,length):
		self._struct = scores
		self._n = length
		self._alphabet = scores.keys()
		self._alphabet.sort()
		#self._scores = {  for kmer in itertools.product( self._alphabet, repeat = self._n ) }	
	
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

	def score(self, seq):
		#assert len(seq) == self.length()
		return sum( self._struct[seq[i]][i] for i in xrange( len( seq ) ) )
		s = 0
		for i in range( len( seq ) ):
			s += self._struct[ seq[i] ][i]
		return s
		#return sum( [ self._struct[seq[i]][i] for i in range(len(seq)) ] )

	def scoreAll(self, s ):
		assert len(s) >= self.length()
		#if tPool != None:
		#	return tPool.map( Score(self,s) , xrange( len(s) - self.length() + 1 ) )
		struct = self._struct
		n = len(s) - len(self) + 1
		m = xrange( len( self ) )
		return ( sum( struct[s[i+j]][j] for j in m ) for i in xrange( n ) )

	def occurs( self, seq, cutoff = 0 ):
		assert len(seq) >= self.length()
		for i in range( len(seq) - self.length() + 1 ):
			if self.score( seq[i:i+self.length()] ) > cutoff:
				return True
		return False

	def occurrences( self, seq, cutoff = 0 ):
		scores = self.scoreAll( seq )
		occ = 0
		for s in scores:
			if s > cutoff: occ += 1
		return occ
	
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

class Score(object):
	def __init__(self, pwm, s):
		self.pwm = pwm
		self.s = s
	def __call__(self, i):
		return self.pwm.score( self.s[i:i+len(self.pwm)] )

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
	denom = float( len(seq) ) + addedCounts
	for (k,n) in counts.iteritems():
		counts[k] = float( n + pseudocount ) / denom
	return counts

class Sample(object):
	def __init__(self, seq, pwm):
		self.seq = seq
		self.pwm = pwm
	def __call__(self, _):
		global S
		global Scoring
		s = time.clock()
		shuf = ''.join( random.sample( self.seq, len(self.seq) ) )
		S += time.clock() - s
		s = time.clock()
		scores = self.pwm.scoreAll( shuf )
		Scoring += time.clock() - s
		return scores

import time
T = 0
S = 0
Sorting = 0
Scoring = 0

def shuffleAndScore( pwm, shuf ):
	#global S
	#global Scoring
	#s2 = time.clock()
	random.shuffle( shuf )
	#S += time.clock() - s2
	#s2 = time.clock()
	shufScores = pwm.scoreAll( shuf )
	#Scoring += time.clock() - s2
	#print [ s for s in shufScores ]
	return shufScores

def nullScoreDistribution( pwm, seq, reps = 1000, tPool = None ):
	#global T
	#global S
	#global Sorting
	#global Scoring
	#s = time.clock()
	#if tPool != None:
	#	scores = tPool.map( Sample( seq, pwm ), xrange( reps ) )
	#	scores = list( itertools.chain.from_iterable( scores ) ) 
	#else:
	shuf = list( seq )
	
	#s2 = time.clock()		
	scores = list( itertools.chain( *(shuffleAndScore( pwm, shuf ) for _ in xrange( reps ) ) ) )
	#print len(scores)
	#Scoring += time.clock() - s2
	
	#s3 = time.clock()
	scores.sort()
	#Sorting += time.clock() - s3
	#T += time.clock() - s
	return scores

def scoreCutoff( nullDist, pval ):
	i = float( len( nullDist ) ) -  pval * len( nullDist )
	upper = int( math.ceil( i ) ) 
	lower = int( math.floor( i ) )
	if upper == lower:
		return nullDist[upper]
	wl = float( upper ) - i
	wu = i - float( lower )
	return wu * nullDist[upper] + wl * nullDist[lower]

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
