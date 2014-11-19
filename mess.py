import sys
import random
import PyPositionWeightMatrix as PWM
import Parser
import Enrichment
import multiprocessing

MOTIF_FLAG = '-m'
SEQ_FLAG = '-s'
CUTOFF_FLAG = '-p'
THREADS_FLAG = '-t'

#returns ( [motif files], [Sequence set files], log likelihood cutoff )
def parseArgs( args ):
	motifFiles = []
	seqFiles = []
	cutoff = 0.05
	pool = None
	cur = None
	i = 0
	while i < len( args ):
		if args[i] == MOTIF_FLAG:
			cur = motifFiles
		elif args[i] == SEQ_FLAG:
			cur = seqFiles
		elif args[i] == CUTOFF_FLAG:
			cur = None
			i += 1
			cutoff = float( args[i] )	
		elif args[i] == THREADS_FLAG:
			cur = None
			i += 1
			pool = multiprocessing.Pool( processes = int( args[i] ) )
		else:
			cur.append( args[i] )
		i += 1
	return motifFiles, seqFiles, cutoff, pool


def unzip( xs ):
	return ( [ a for a,b in xs ] , [ b for a,b in xs ] )

class Label(object):

	def __init__(self, seqSet, bgs, pwm, cutoff):
		self.seqSet = seqSet
		self.bgs = bgs
		#print pwm
		self.pwm = pwm
		#print self.pwmI
		self.cutoff = cutoff

	def __call__(self, i):
		#print self.pwm
		(fc, s) = self.seqSet[i]
		#print self.pwm
		label = labelSequence( s, self.bgs[i], self.pwm, self.cutoff )
		#print >> sys.stderr, PWM.T,PWM.S,PWM.Sorting,PWM.Scoring
		return (fc,label)

def labelSequence( s, bg, pwm, cutoff ):
	#print pwm
	#adjust the PWM to the sequence background nucleotide frequency
	pwmAdj = pwm.adjustToBG( bg )
	#print pwmAdj
	#estimate log likelihood cutoff from sampled null distribution
	nullDist = PWM.nullScoreDistribution( pwmAdj, s )
	scoreCutoff = PWM.scoreCutoff( nullDist, cutoff )
	
	label = 1 if pwmAdj.occurs( s, scoreCutoff ) else -1
	return label

def motifEnrichment( seqSet, pwm, cutoff, bgs = None, isSorted = False, tPool = None ):
	
	#print pwm

	if bgs == None:
		bgs = [ PWM.nucleotideFrequency( s ) for (_,s) in seqSet ]
	assert len( seqSet ) == len( bgs )
	
	if tPool != None:
		#args = [ ( seqSet[i][1], bgs[i], pwm, cutoff ) for i in xrange( len( seqSet ) ) ]
		#data = tPool.map( labelSequence, args )
		#data = [ ( seqSet[i][0], data[i] ) for i in xrange( len( data ) ) ]
		data = tPool.map( Label( seqSet, bgs, pwm, cutoff ), xrange( len( seqSet ) ) )
	else:
		data = []
		for i in range( len( seqSet ) ):
			(fc, s) = seqSet[i]
			label = labelSequence( s, bgs[i], pwm, cutoff )
			#fc,label = Label( seqSet, bgs, pwm, cutoff ) ( i )
			data.append( (fc, label) )
	
	if not isSorted:
		data.sort( reverse = True )
	
	( w, l ) = unzip( data )
	#( sig, es, posScores ) = Enrichment.enrichment( w, l, preSorted = True )	
	return Enrichment.enrichment( w, l, preSorted = True, tPool = tPool )

def translate( seq, alphabet ):
	return [ alphabet[c] for c in seq ]

def parseSeqFile( seqFile, alphabet ):
	seqs = [ (foldChange, translate(s,alphabet)) for (_,foldChange,s) in Parser.parse( open( seqFile, 'r' ) ) ] 
	seqs.sort( reverse = True )
	
	return seqs
		
class EnrichmentCalculator(object):

	def __init__(self, seqSet, cutoff, bgs ):
		self.seqSet = seqSet
		self.cutoff = cutoff
		self.bgs = bgs

	def __call__(self, pwm):
		return motifEnrichment( self.seqSet, pwm, self.cutoff, bgs = self.bgs, isSorted = True )

def nucleotideFrequency( seqEntry ):
	_,s = seqEntry
	return PWM.nucleotideFrequency(s)

def alphabetListToDict( alph ):
	return { c : i for ( i, c ) in enumerate( alph ) }

def computeMotifEnrichments( seqFiles, pwms, cutoff, alphabet, tPool = None ):
	alphDict = alphabetListToDict( alphabet )
	for f in seqFiles:
		seqSet = parseSeqFile( f, alphDict )
		if tPool != None:
			bgs = tPool.map( nucleotideFrequency, seqSet )
		else:
			bgs = [ PWM.nucleotideFrequency( s ) for (_,s) in seqSet ]
		
		report =  [ motifEnrichment( seqSet, pwm, cutoff, bgs = bgs, isSorted = True, tPool = tPool ) for pwm in pwms  ]
		yield f , report

def printReport( motifFiles, report, header = True ):
	assert len( motifFiles ) == len( report )
	assert len( report ) > 0
	
	#print header
	if header:
		print 'PWM', 'Significance', 'ES',
		( _, _, pos ) = report[0]
		for i in range( len( pos ) ):
			print i+1,
		print ''
	
	#print report
	for i in range( len( report ) ):
		( sig, es, pos ) = report[i]
		print motifFiles[i], sig, es,
		for s in pos:
			print s,
		print ''
	

def main( args ):
	#do stuff
	motifFiles, seqFiles, cutoff, pool = parseArgs( args )
	pwms = [ PWM.parsePFM( open( f, 'r' ) ) for f in motifFiles ]
	#for pwm in pwms:
	#	print pwm
	alphabet = pwms[0].alphabet()
	
	#pool = multiprocessing.pool.ThreadPool( processes = 4 )
	#pool = multiprocessing.Pool( processes = 4 )

	for f, report in computeMotifEnrichments( seqFiles, pwms, cutoff, alphabet, pool ):
		print f
		printReport( motifFiles, report )	

	#print >> sys.stderr, 'Null sequence dist time = ',PWM.T,PWM.S,PWM.Sorting,PWM.Scoring

	"""
	import time
	t = time.clock()
	for f, report in computeMotifEnrichments( seqFiles, pwms, cutoff ):
		print f
		printReport( motifFiles, report )	
	t = time.clock() - t
	print 'Unthreaded time taken = ',t
	
	t = time.clock()
	for f, report in computeMotifEnrichments( seqFiles, pwms, cutoff, pool ):
		print f
		printReport( motifFiles, report )	
	t = time.clock() - t
	print 'Threaded time taken = ',t
	"""

	"""
	testSeq = 'TCGAATATTACGA'
	for p in pwms:
		print p
		print p.scoreAll( testSeq )
		print occurrences( p, testSeq, 0 )
	"""	
	
	"""
	random.seed(2)
	
	testWeights = [ 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5 ]
	testLabels = [ 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, -1 ]
	print Enrichment.enrichment( testWeights, testLabels, rho = 0 )
	print Enrichment.enrichment( testWeights, testLabels, rho = 1 )
	
	testLabels = [ -1 * x for x in testLabels ]
	print Enrichment.enrichment( testWeights, testLabels , rho = 0 )
	print Enrichment.enrichment( testWeights, testLabels, rho = 1 )
	
	
	testWeights = [ 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5 ]
	testLabels = [ 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1 ]
	print Enrichment.enrichment( testWeights, testLabels , rho = 0 )
	print Enrichment.enrichment( testWeights, testLabels, rho = 1 )
	
	testLabels = [ -1 * x for x in testLabels ]
	print Enrichment.enrichment( testWeights, testLabels , rho = 0 )
	print Enrichment.enrichment( testWeights, testLabels, rho = 1 )
	
	testWeights = [ 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5 ]
	testLabels = [ -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1 ]
	print Enrichment.enrichment( testWeights, testLabels , rho = 0 )
	print Enrichment.enrichment( testWeights, testLabels, rho = 1 )
	
	testLabels = [ -1 * x for x in testLabels ]
	print Enrichment.enrichment( testWeights, testLabels , rho = 0 )
	print Enrichment.enrichment( testWeights, testLabels, rho = 1 )
	"""

if __name__ == '__main__':
	main( sys.argv[1:] )
