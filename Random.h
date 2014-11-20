/*

Author: Tristan Bepler

*/


#ifndef RANDOM_H
#define RANDOM_H

#include <stdlib.h>
#include <stdio.h>
#include "Arrays.h"

//generates a random double between 0 and 1
inline static double randomDouble( void )
{
	return (double)random( )/(double)RAND_MAX;
}

//courtesy of http://stackoverflow.com/questions/2509679/how-to-generate-a-random-number-from-within-a-range
// Assumes 0 <= range <= RAND_MAX
// Returns in the half-open interval [0, max]
inline static long randomRange(long min, long max) {
	
	max = max - min;

	unsigned long
	// max <= RAND_MAX < ULONG_MAX, so this is okay.
	numBins = (unsigned long) max + 1,
    	numRand = (unsigned long) RAND_MAX + 1,
    	binSize = numRand / numBins,
    	defect   = numRand % numBins;

  	long x;
  	// This is carefully written not to overflow
  	while (numRand - defect <= (unsigned long)(x = random()));

  	// Truncated division is intentional
  	return min + x/binSize;
}

//shuffles the array in place
inline static void shuffle( void* array, size_t len, size_t size )
{
	char* arr = (char*) array;
	size_t i;
	size_t j;
	size_t n = len - 1;
	for( i = 0 ; i < n ; ++i )
	{
		j = randomRange( i, n );
		//printf( "%d\n", j );
		swap( arr + i*size, arr + j*size, size );
	}
}




#endif /* RANDOM_H */
