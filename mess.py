import sys
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

def occurrences( pwm, seq, cutoff ):
	scores = pwm.scoreAll( seq )
	#print scores
	occ = 0
	for s in scores:
		if s > cutoff: occ += 1
	return occ

def enrichment( pwm, seqFile, logLikelihoodCutoff ):
	seqs = [ (foldChange, s) for (_,foldChange,s) in Parser.parse( open( seqFile, 'r' ) ) ] 
	seqs.sort( revers = True )

	motif = {}
	for (foldChange, seq) in seqs:
		bg = PWM.nucleotideFrequency( seq )
		pwmAdj = pwm.adjustToBG( bg )
		motif[ seq ] = occurrences( pwmAdj, seq, logLikelihoodCutoff ) > 0

	
		

def main( args ):
	#do stuff
	motifFiles, seqFiles, cutoff = parseArgs( args )
	pwms = [ PWM.parsePFM( open( f, 'r' ) ) for f in motifFiles ]
	
	"""
	testSeq = 'TCGAATATTACGA'
	for p in pwms:
		print p
		print p.scoreAll( testSeq )
		print occurrences( p, testSeq, 0 )
	"""	

	testWeights = [ 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5 ]
	testLabels = [ 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, -1 ]
	print Enrichment.enrichment( testWeights, testLabels )
	
	testWeights = [ 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5 ]
	testLabels = [ 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1 ]
	print Enrichment.enrichment( testWeights, testLabels )


if __name__ == '__main__':
	main( sys.argv[1:] )
