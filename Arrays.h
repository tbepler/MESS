/*

Author: Tristan Bepler

*/

#ifndef ARRAYS_H
#define ARRAYS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

inline static void swap( void* i, void* j, size_t size )
{
	char temp[size];
	memcpy(temp, i, size);
	memcpy(i, j, size);
	memcpy(j, temp, size);
}

inline static size_t partition( void* array, size_t size, size_t low, size_t up, int (*compare)( const void*, const void* ) )
{
	char* arr;
	arr = (char*) array;
	char* pivot = arr + ( low + ((( up - low ) / size) / 2 ) * size );
	char x[size];
	memcpy( x, pivot, size );

	//swap pivot to the end
	swap( pivot, arr + up, size );	

	size_t i = low;
	size_t j;
	for( j = low ; j < up ; j += size )
	{
		//printf( "j = %d\n", j );
		if( compare( arr + j, &x ) <= 0 )
		{
			swap( arr + i , arr + j, size );
			i += size;
		}
	}
	//swap pivot into correct position
	swap( arr + i, arr + up, size );
	return i;
}

inline static void printDoubleArray( double* arr, size_t len )
{
	int i;
	for( i = 0 ; i < len ; ++i )
	{
		printf( "%f", arr[i] );
		if( i != len -1 )
		{
			printf( ", " );
		}
	}
	printf( "\n" );
}	

inline static void quicksort( void* array, size_t len, size_t size, int (*compare)( const void*, const void* ) )
{
	
	char* arr = (char*) array;	
	
	size_t low = 0;
	size_t up = (len - 1) * size;

	size_t stack[ len - 1 ];
	ssize_t top = 0;
	
	stack[ top ] = low;
	stack[ ++top ] = up;

	size_t part;
	while( top >= 0 )
	{
		up = stack[ top-- ];
		low = stack[ top-- ];

		part = partition( arr, size, low, up, compare );

		if( part > low + size )
		{
			stack[ ++top ] = low;
			stack[ ++top ] = part - size;
		}
		if( part + size < up )
		{
			stack[ ++top ] = part + size;
			stack[ ++top ] = up;
		}
	}
}

inline static int isSorted( const void* array, size_t len, size_t size, int (*compare)( const void*, const void* ) )
{
	char* arr = (char*) array;
	size_t i;
	size_t end = ( len - 1 ) * size;
	for( i = 0 ; i < end ; i += size )
	{
		if( compare( arr + i, arr + i + size ) > 0 )
		{
			return 0;
		}
	}
	return 1;
}

inline static int compareDoubleAscending( const void* i, const void* j )
{
	double a = * (double*) i;
	double b = * (double*) j;
	if( a < b )
	{
		return -1;
	}
	if( a > b )
	{
		return 1;
	}
	return 0;
}

inline static int compareDoubleDescending( const void* i, const void* j )
{
	double a = * (double*) i;
	double b = * (double*) j;
	if( a > b )
	{
		return -1;
	}
	if( a < b )
	{
		return 1;
	}
	return 0;
}





#endif /* ARRAYS_H */

