import random
import copy
import bisect

""" returns significance, enrichment score, rank scores
labels must be 1 and -1 where 1 is in S and -1 is not
in S """
def enrichment( weights, labels, null = None ):
	( es, scores ) = enrichmentScore( weights, labels )
	if null == None:
		null = nullDistribution( weights, labels )
	s = significance( es, null )
	return ( s, es, scores )

def significance( es, null ):
	i = bisect.bisect_left( null, es )
	return 1.0 - float(i) / float( len( null ) )

"""labels must be either 1 or -1 for in the set
or not in the set """
def enrichmentScore( weights, labels ):
	assert len(weights) == len(labels)
	data = zip( weights, labels )
	data.sort( reverse = True )
	
	scores = []
	cur = 0
	for i in range( len( data ) ):
		( w, label ) = data[i]
		cur += w * label
		scores.append( cur )

	return max(scores), scores

def nullDistribution( weights, labels, samples = 1000 ):
	scores = []
	weights = sorted( weights, reverse = True )
	labels = copy.copy( labels )
	for _ in range( samples ):
		random.shuffle(labels)
		( es, _ ) = enrichmentScore( weights, labels )
		scores.append( es )
	scores.sort()
	return scores


