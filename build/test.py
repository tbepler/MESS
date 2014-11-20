import PositionWeightMatrix as PWM
import math
import copy

def transpose( x ):
	return map( list, zip(*x) )

M = transpose( [
	[ 0.1, 0.25, 0.25, 0.8, 0.9, 0.1 ],
	[ 0.2, 0.1, 0.15, 0.05, 0.01, 0.5 ],
	[ 0.6, 0.35, 0.2, 0.05, 0.04, 0.2 ],
	[ 0.1, 0.3, 0.4, 0.1, 0.05, 0.2 ]
	] )

S = [ 0, 0, 1, 2, 0, 3 ]
S2 = [ 1, 2, 1, 0, 3, 1, 2, 0, 0, 0, 1, 2, 0, 3 ]

def log( x ):
	if x == 0: return float('-inf')
	return math.log( x )

def main():
	#testing testing 1 2 3
	global M
	M = [ [ log( x ) for x in xs] for xs in M ]
	pwm = PWM.PositionWeightMatrix( M )
	print pwm.score(S)
	print pwm.scoreAll( S2 )
	print len( pwm )
	print pwm.length()
	for row in pwm:
		print row
	print ''
	c = copy.copy( pwm )	
	for row in c:
		print row
	print ''
	c = copy.deepcopy( pwm )
	for row in c:
		print row
	print ''
	print pwm.nullScoreDistribution( S2, 10 )
	
if __name__ == '__main__':
	main()

