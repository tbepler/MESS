import random
import copy
import bisect
import multiprocessing

""" returns significance, enrichment score, rank scores
labels must be 1 and -1 where 1 is in S and -1 is not
in S
"""
def enrichment( weights, labels, rho = 1, nullDist = None, samples = 1000, preSorted = False, tPool = None ):
	( es, scores ) = enrichmentScore( weights, labels, rho, preSorted )
	if nullDist == None:
		nullDist = nullDistribution( weights, labels, samples, rho, preSorted, tPool )
	s = significance( es, nullDist )
	return ( s, es, scores )

def significance( es, null ):
	if es < 0:
		i = bisect.bisect( null, es )
		return float(i) / float( len( null ) )
	else:
		i = bisect.bisect_left( null, es )
		return 1.0 - float(i) / float( len( null ) )

"""labels must be either 1 or -1 for in the set
or not in the set
rho controls step weight as w^rho"""
def enrichmentScore( weights, labels, rho = 1, preSorted = False ):
	assert len(weights) == len(labels)
	data = zip( weights, labels )
	if not preSorted:
		data.sort( reverse = True )
	
	expected = float( sum( labels ) ) / float( len (labels ) )
	scores = []
	cur = 0
	for i in range( len( data ) ):
		( w, label ) = data[i]
		cur += ( label - expected ) * w ** rho
		scores.append( cur )

	minS = min( scores )
	maxS = max( scores )
	
	es = minS if abs(minS) > maxS else maxS

	return es, scores

class Score(object):
	def __init__(self, weights, labels, rho, preSorted):
		self.w = weights
		self.l = labels
		self.r = rho
		self.s = preSorted
	def __call__(self, _):
		labels = random.sample( self.l, len( self.l ) )
		(es, _) = enrichmentScore( self.w, labels, self.r, self.s)
		return es

def nullDistribution( weights, labels, samples = 1000, rho = 1, preSorted = False, tPool = None ):
	scores = []
	if not preSorted:
		weights = sorted( weights, reverse = True )
	labels = copy.copy( labels )
	if tPool != None:
		scores = tPool.map( Score( weights, labels, rho, preSorted ), xrange( samples ) )
	else:
		for _ in range( samples ):
			random.shuffle(labels)
			( es, _ ) = enrichmentScore( weights, labels, rho, preSorted = True )
			scores.append( es )
	scores.sort()
	return scores


