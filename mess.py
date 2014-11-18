import sys
import random
import PositionWeightMatrix as PWM
import Parser
import Enrichment

MOTIF_FLAG = '-m'
SEQ_FLAG = '-s'
LOG_LIKELIHOOD_FLAG = '-l'

#returns ( [motif files], [Sequence set files], log likelihood cutoff )
def parseArgs( args ):
	motifFiles = []
	seqFiles = []
	cutoff = 0
	cur = None
	i = 0
	while i < len( args ):
		if args[i] == MOTIF_FLAG:
			cur = motifFiles
		elif args[i] == SEQ_FLAG:
			cur = seqFiles
		elif args[i] == LOG_LIKELIHOOD_FLAG:
			cur = None
			i += 1
			cutoff = float( args[i] )	
		else:
			cur.append( args[i] )
		i += 1
	return motifFiles, seqFiles, cutoff


def unzip( xs ):
	return ( [ a for a,b in xs ] , [ b for a,b in xs ] )

def motifEnrichment( seqSet, pwm, logLikelihoodCutoff, bgs = None, isSorted = False ):
	if bgs == None:
		bgs = [ PWM.nucleotideFrequency( s ) for (_,s) in seqSet ]
	assert len( seqSet ) == len( bgs )
	
	data = []
	for i in range( len( seqSet ) ):
		(fc, s) = seqSet[i]
		pwmAdj = pwm.adjustToBG( bgs[i] )
		#print pwmAdj.scoreAll( s )
		label = 1 if pwmAdj.occurs( s, logLikelihoodCutoff ) else -1
		#print label
		data.append( (fc, label) )
	
	if not isSorted:
		data.sort( reverse = True )
	
	( w, l ) = unzip( data )
	#( sig, es, posScores ) = Enrichment.enrichment( w, l, preSorted = True )	
	return Enrichment.enrichment( w, l, preSorted = True )

def parseSeqFile( seqFile ):
	seqs = [ (foldChange, s) for (_,foldChange,s) in Parser.parse( open( seqFile, 'r' ) ) ] 
	seqs.sort( reverse = True )
	return seqs
		
def computeMotifEnrichments( seqFiles, pwms, logLikelihoodCutoff ):
	reports = []
	for f in seqFiles:
		seqSet = parseSeqFile( f )
		bgs = [ PWM.nucleotideFrequency( s ) for (_,s) in seqSet ]
		es = [ motifEnrichment( seqSet, pwm, logLikelihoodCutoff, bgs, True ) for pwm in pwms  ]
		reports.append( es )
	return reports	

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
	motifFiles, seqFiles, cutoff = parseArgs( args )
	pwms = [ PWM.parsePFM( open( f, 'r' ) ) for f in motifFiles ]
	
	reports = computeMotifEnrichments( seqFiles, pwms, cutoff )
	for i in range( len( reports ) ):
		print seqFiles[i]
		printReport( motifFiles, reports[i] )	
	
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
