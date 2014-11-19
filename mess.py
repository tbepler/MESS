import sys
import random
import PyPositionWeightMatrix as PWM
import Parser
import Enrichment
import multiprocessing
import os

MOTIF_FLAG = '-m'
SEQ_FLAG = '-s'
CUTOFF_FLAG = '-p'
THREADS_FLAG = '-t'
DETAILS_FLAG = '-d'

#returns ( [motif files], [Sequence set files], log likelihood cutoff, thread pool, details dir )
def parseArgs( args ):
	motifFiles = []
	seqFiles = []
	cutoff = 0.05
	pool = None
	detailDir = None
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
		elif args[i] == DETAILS_FLAG:
			cur = None
			i += 1
			detailDir = args[i]
		else:
			cur.append( args[i] )
		i += 1
	return motifFiles, seqFiles, cutoff, pool, detailDir


def unzip( xs ):
	if( len( xs ) < 1 ): return xs
	n = len( xs[0] )
	if n == 1:
		return xs
	elif n == 2:
		return unzip2( xs )
	elif n == 3:	
		return unzip3( xs )
	elif n == 4:
		return unzip4( xs )
	elif n == 5:
		return unzip5( xs )
	raise InputError( xs, "Invalid unzip tuple size, " + str(n) ) 

def unzip2( xs ):
	return ( [ a for a,b in xs ] , [ b for a,b in xs ] )

def unzip3( xs ):
	return ( [ a for a,b,c in xs ] , [ b for a,b,c in xs ], [ c for a,b,c in xs ] )

def unzip4( xs ):
	return ( 
		[ a for a,b,c,d in xs ] ,
		[ b for a,b,c,d in xs ] ,
		[ c for a,b,c,d in xs ] ,
		[ d for a,b,c,d in xs ]
		)

def unzip5( xs ):
	return ( 
		[ a for a,b,c,d,e in xs ] ,
		[ b for a,b,c,d,e in xs ] ,
		[ c for a,b,c,d,e in xs ] ,
		[ d for a,b,c,d,e in xs ] ,
		[ e for a,b,c,d,e in xs ]
		)

class Label(object):

	def __init__(self, seqSet, bgs, pwm, cutoff):
		self.seqSet = seqSet
		self.bgs = bgs
		self.pwm = pwm
		self.cutoff = cutoff

	def __call__(self, i):
		(fc, n, s) = self.seqSet[i]
		(label, scoreCutoff, count) = labelSequence( s, self.bgs[i], self.pwm, self.cutoff )
		return (fc, n, label, scoreCutoff, count)

def labelSequence( s, bg, pwm, cutoff ):
	#adjust the PWM to the sequence background nucleotide frequency
	pwmAdj = pwm.adjustToBG( bg )
	#estimate log likelihood cutoff from sampled null distribution
	nullDist = PWM.nullScoreDistribution( pwmAdj, s )
	scoreCutoff = PWM.scoreCutoff( nullDist, cutoff )
	count = pwmAdj.occurrences( s, scoreCutoff )

	label = 1 if count > 0 else -1
	return label, scoreCutoff, count 

def motifEnrichment( seqSet, pwm, cutoff, bgs = None, isSorted = False, tPool = None ):
	
	#print pwm

	if bgs == None:
		bgs = [ PWM.nucleotideFrequency( s ) for (_,s) in seqSet ]
	assert len( seqSet ) == len( bgs )
	
	if tPool != None:
		data = tPool.map( Label( seqSet, bgs, pwm, cutoff ), xrange( len( seqSet ) ) )
	else:
		data = []
		for i in range( len( seqSet ) ):
			(fc, n, s) = seqSet[i]
			(label,scoreCutoff,count) = labelSequence( s, bgs[i], pwm, cutoff )
			data.append( (fc, n, label, scoreCutoff, count) )
	
	if not isSorted:
		data.sort( reverse = True )
	
	( weights, names, labels, scoreCutoffs, counts ) = unzip( data )
	( sig, es, posScores ) = Enrichment.enrichment( weights, labels, preSorted = True, tPool = tPool )	
	return  sig, es, posScores, names, labels, scoreCutoffs, counts

def translate( seq, alphabet ):
	return [ alphabet[c] for c in seq ]

def parseSeqFile( seqFile, alphabet ):
	f = open( seqFile, 'r' )
	seqs = [ (foldChange, name, translate(s,alphabet)) for (name,foldChange,s) in Parser.parse( f ) ]
	f.close() 
	seqs.sort( reverse = True )
	
	return seqs

def nucleotideFrequency( seqEntry ):
	_,_,s = seqEntry
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

def toDetailFile( detailDir, motifFileName, seqFile = None ):
	mName = os.path.splitext( os.path.basename( motifFileName ) )[0]
	if seqFile != None:
		sName = os.path.splitext( os.path.basename( seqFile ) )[0]
		mName = sName + '_' + mName
	fName = os.path.join( detailDir, mName ) + '_details.txt'
	if not os.path.exists( detailDir ):
		os.makedirs( detailDir )
	return open( fName, 'w' )
	

def printReport( motifFiles, report, header = True, detailDir = None, seqFile = None ):
	assert len( motifFiles ) == len( report )
	assert len( report ) > 0
	#report is ( sig, es, posScores, names, labels, scoreCutoffs, counts )

	#print header
	if header:
		print 'PWM', 'Significance', 'ES',
		( _, _, pos, _, _, _, _ ) = report[0]
		for i in range( len( pos ) ):
			print i+1,
		print ''
	
	#print report
	for i in range( len( report ) ):
		( sig, es, pos, _, _, _, _ ) = report[i]
		print motifFiles[i], sig, es,
		for s in pos:
			print s,
		print ''

	#write details to file if specified
	if detailDir != None:
		for i in range( len( report ) ):
			f = toDetailFile( detailDir, motifFiles[i], seqFile )
			f.write('Name ScoreCutoff Occurrences Label')
			( _, _, _, names, labels, scoreCutoffs, counts ) = report[i]
			for j in range( len( names ) ):
				f.write( '\n' + ' '.join( map( str, [ names[j], scoreCutoffs[j], counts[j], labels[j] ] ) ) )
			f.close()

def main( args ):
	#do stuff
	motifFiles, seqFiles, cutoff, pool, detailDir = parseArgs( args )
	pwms = [ PWM.parsePFM( open( f, 'r' ) ) for f in motifFiles ]
	for pwm in pwms:
		print pwm
	alphabet = pwms[0].alphabet()
	
	#pool = multiprocessing.pool.ThreadPool( processes = 4 )
	#pool = multiprocessing.Pool( processes = 4 )

	for f, report in computeMotifEnrichments( seqFiles, pwms, cutoff, alphabet, pool ):
		print f
		printReport( motifFiles, report, detailDir = detailDir, seqFile = f )	

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
