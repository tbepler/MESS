import sys
import PositionWeightMatrix as PWM

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

def occurrences( pwm, seq, logLikelihoodCutoff ):
	bg = PWM.nucleotideFrequency( seq )
	pwmAdj = pwm.adjustToBG( bg )
	scores = pwmAdj.scoreAll( seq )
	#print scores
	occ = 0
	for s in scores:
		if s > logLikelihoodCutoff: occ += 1
	return occ


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

if __name__ == '__main__':
	main( sys.argv[1:] )
