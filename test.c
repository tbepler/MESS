#include <stdlib.h>
#include <stdio.h>
#include "Arrays.h"
#include "Random.h"

static const int LEN = 1000000;
//static double array[LEN];

int main( int argc, const char* argv[] )
{
	//double array[] = { 1, 1.5, -53, 298,-139, -13841, -0.2882, 1.814, 49292.131901, -1814.8148, 0, 0.131 };
	//size_t LEN = 12;
	
	//const int REAL_BIG = 10000000;
	printf( "Allocating array of size %d...\n",LEN );
	double* array = malloc( LEN * sizeof(double) );
	//size_t LEN = REAL_BIG;
	
	const double mag = 100000;
	const double adj = mag / 2;	
	
	printf( "Initializing array...\n" );
	int i;
	for( i = 0 ; i < LEN ; ++i )
	{
		array[i] = randomDouble() * mag - adj;
	}

	//printDoubleArray( array, LEN );
	printf( "Is sorted? %d\n", isSorted( array, LEN, sizeof(double), &compareDoubleAscending ) );
	printf( "Sorting...\n");
	quicksort( array, LEN, sizeof(double), &compareDoubleAscending );
	printf( "Is sorted? %d\n", isSorted( array, LEN, sizeof(double), &compareDoubleAscending ) );
	//printDoubleArray( array, LEN );
	printf( "Shuffling...\n");
	shuffle( array, LEN, sizeof(double) );
	printf( "Is sorted? %d\n", isSorted( array, LEN, sizeof(double), &compareDoubleAscending ) );
	//printDoubleArray( array, LEN );
	printf( "Sorting again...\n" );	
	quicksort( array, LEN, sizeof(double), &compareDoubleDescending );
	printf( "Is sorted? %d\n", isSorted( array, LEN, sizeof(double), &compareDoubleDescending ) );
	//printDoubleArray( array, LEN );
	
	
}


