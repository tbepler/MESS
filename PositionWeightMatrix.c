/*

Author: Tristan Bepler


*/

#include <Python.h>
#include <structmember.h>
#include "Arrays.h"
#include "Random.h"

typedef struct {
	PyObject_HEAD
	
	double** matrix;
	size_t n;
	size_t m;

} PositionWeightMatrix;

static PyObject* PositionWeightMatrix_new( PyTypeObject* type, PyObject* args, PyObject* kwds )
{
	PositionWeightMatrix* self;
	self = (PositionWeightMatrix*) type -> tp_alloc( type, 0 );
	self -> matrix = NULL;
	self -> n = self -> m = 0;
	return ( PyObject * ) self;
}

inline static void freeMatrix( double** matrix, int n )
{
	if( matrix )
	{
		//free each row of the matrix
		int i;
		for( i = 0 ; i < n ; ++i )
		{
			free( matrix[i] );
		}
		//free the top level pointer
		free( matrix );
	}
}

static void PositionWeightMatrix_dealloc( PositionWeightMatrix* self )
{
	int n = self -> n;
	double** matrix = self -> matrix;
	//free matrix
	freeMatrix( matrix, n );
	//free python object
	self -> ob_type -> tp_free( (PyObject*) self );
}

inline static void resetPositionWeightMatrix( PositionWeightMatrix* self )
{
	//first free matrix
	freeMatrix( self->matrix, self->n );
	//reset contents
	self->matrix = NULL;
	self->n = 0;
	self->m = 0;
}

inline static int initializeMatrix( PositionWeightMatrix* self, PyObject* list )
{
	//get the number of elements in the list
	int n = PyList_Size( list );
	if( n < 0 ) return -1; //should raise error

	if( n == 0 )
	{
		resetPositionWeightMatrix( self );
		return 0; //no need to proceed
	}

	//precheck rows before allocation memory
	int m = -1;
	PyObject* row;
	int i;
	for( i = 0 ; i < n ; ++i )
	{
		row = PyList_GetItem( list, i );
		if( PyList_Check( row ) )
		{
			if( m < 0 )
			{
				//assign row length
				m = PyList_Size( row );
			}
			else if( m != PyList_Size( row ) ) //ensure all rows are same length
			{
				//TODO raise error here
				return -1;
			}

		}
		else
		{
			//TODO raise error here
			return -1;
		}
	}
	if( m == 0 )
	{
		resetPositionWeightMatrix( self );
		return 0; //no need to proceed
	}

	//allocate matrix with n rows and m columns
	double** matrix = malloc( sizeof( double* ) * n );
	for( i = 0 ; i < n ; ++i )
	{
		matrix[i] = malloc( sizeof( double ) * m );
	}

	//fill matrix with values from the list rows
	double d;
	int j;
	for( i = 0 ; i < n ; ++i )
	{
		row = PyList_GetItem( list, i );
		for( j = 0 ; j < m ; ++j )
		{
			d = PyFloat_AsDouble( PyList_GetItem( row, j ) );
			//check for errors in converting to float
			if( PyErr_Occurred() )
			{
				//free memory and return
				freeMatrix( matrix, n );
				return -1;
			}
			//assign matrix entry
			matrix[i][j] = d;	
		}
	}

	//assign matrix, n, and m into self
	//first free current matrix in self - just in case
	freeMatrix( self->matrix, self->n );
	self->matrix = matrix;
	self->n = n;
	self->m = m;

	//return success
	return 0;
	
}

static int PositionWeightMatrix_init( PositionWeightMatrix* self, PyObject* args, PyObject* kwds )
{
	PyObject* listObj;

	static char *kwlist[] = { "matrix", NULL };
	
	if( !PyArg_ParseTupleAndKeywords( args, kwds, "O!", kwlist, &PyList_Type, &listObj ) )
		return -1;
	
	if( listObj )
	{
		return initializeMatrix( self, listObj );
	}
	
	//arg is not list, so error
	return -1;
	
}

inline static PyObject* parseTupleToFastSequence( PyObject* args )
{
	PyObject* seq;
	if( !PyArg_ParseTuple( args, "O", &seq ) ) return NULL;
	return PySequence_Fast( seq, "argument must be iterable" );
}

inline static double score( double** matrix, int n, int m, PyObject** elems )
{
	double score = 0;
	int i;
	int j;
	for( i = 0 ; i < n ; ++i )
	{
		j = ( int ) PyInt_AsLong( elems[ i ] );
		//should check for error here
		//if( PyErr_Occurred() )
		score += matrix[i][j];
	}

	return score;
}

inline static double scoreFromIndex( double** matrix, int n, int m, PyObject** elems, int start )
{
	double score = 0;
	int i;
	int j;
	for( i = 0 ; i < n ; ++i )
	{
		j = ( int ) PyInt_AsLong( elems[ i + start ] );
		//should check for error here
		//if( PyErr_Occurred() )
		score += matrix[i][j];
	}

	return score;
}

inline static double scoreLongFromIndex( double** matrix, size_t n, size_t m, long* seq, size_t start )
{
	double score = 0;
	size_t i;
	long j;
	for( i = 0 ; i < n ; ++i )
	{
		j = seq[ i + start ];
		score += matrix[i][j];
	}
	return score;
}

inline static void scoreAllLong( double** matrix, size_t n, size_t m,
				long* seq, size_t len,
				double* scores, size_t start )
{
	size_t end = len - n + 1;
	size_t i;
	for( i = 0 ; i < end ; ++i )
	{
		scores[ start + i ] = scoreLongFromIndex( matrix, n, m, seq, i );
	}
}

static PyObject * PositionWeightMatrix_score( PositionWeightMatrix* self, PyObject* args )
{
	PyObject* seq = parseTupleToFastSequence( args );
	if( !seq ) return NULL;
	
	double** matrix = self -> matrix;
	int n = self -> n;
	int m = self -> m;
	int len = PySequence_Fast_GET_SIZE( seq );
	//check that length is correct
	if( n != len )
	{
		//TODO raise error
		return NULL;
	}
	
	PyObject** elems = PySequence_Fast_ITEMS( seq );
	
	double s = score( matrix, n, m, elems ); 

	return PyFloat_FromDouble(s);
	
}

inline static PyObject* toPyList( double* arr, int n )
{
	PyObject* list = PyList_New( n );
	int i;
	for( i = 0 ; i < n ; ++i )
	{
		PyList_SetItem( list, i, PyFloat_FromDouble( arr[i] ) );
	}
	return list;
}

static PyObject* PositionWeightMatrix_scoreAll( PositionWeightMatrix* self, PyObject* args )
{
	PyObject* seq = parseTupleToFastSequence( args );
	if( !seq ) return NULL;
	
	double** matrix = self -> matrix;
	int n = self -> n;
	int m = self -> m;
	int len = PySequence_Fast_GET_SIZE( seq );
	//check that length is correct
	if( n > len )
	{
		//TODO raise error
		return NULL;
	}
	
	PyObject** elems = PySequence_Fast_ITEMS( seq );
	
	//create new python list and fill it with the scores at each position
	int end = len - n + 1;
	PyObject* list = PyList_New( end );
	int i;
	for( i = 0 ; i < end ; ++i )
	{
		PyList_SetItem( list, i, PyFloat_FromDouble( scoreFromIndex( matrix, n, m, elems, i ) ) );
	}
	return list;
	
}

static PyObject* PositionWeightMatrix_length( PositionWeightMatrix* self )
{
	return PyInt_FromLong( self -> n );
}

static Py_ssize_t PositionWeightMatrix_len( PositionWeightMatrix* self )
{
	return self -> n;
}

static PyObject* PositionWeightMatrix_getItem( PositionWeightMatrix* self, Py_ssize_t i )
{
	//return the jth column of the matrix as a python list
	double** matrix = self -> matrix;
	int n = self -> n;
	int m = self -> m;
	
	if( i < 0 || i >= n )
	{
		//TODO raise error
		return NULL;
	}

	PyObject* list = PyList_New( m );
	int j;
	for( j = 0 ; j < m ; ++j )
	{
		PyList_SetItem( list, j, PyFloat_FromDouble( matrix[i][j] ) );
	}
	return list;
}

static PyObject* PositionWeightMatrix_toList( PositionWeightMatrix* self )
{
	double** matrix = self -> matrix;
	int n = self -> n;
	int m = self -> m;
	
	PyObject* list = PyList_New( n );
	PyObject* row;
	int i;
	int j;
	for( i = 0 ; i < n ; ++i )
	{
		row = PyList_New( m );
		for( j = 0 ; j < m ; ++j )
		{
			PyList_SetItem( row, j, PyFloat_FromDouble( matrix[i][j] ) );
		}
		PyList_SetItem( list, i, row );
	}	
	return list;
}

static PyTypeObject PositionWeightMatrixType;

static PyObject* PositionWeightMatrix_copy( PositionWeightMatrix* self )
{
	PyObject* args = PositionWeightMatrix_toList( self );
	PyObject* tuple = PyTuple_Pack( 1, args );
	PyObject* copy = PyObject_CallObject( (PyObject*) &PositionWeightMatrixType, tuple );
	Py_DECREF( args );
	Py_DECREF( tuple );
	return copy;
	//return (PyObject*) self;
}

static PyObject* PositionWeightMatrix_deepcopy( PositionWeightMatrix* self, PyObject* kywds )
{
	PyObject* args = PositionWeightMatrix_toList( self );
	PyObject* tuple = PyTuple_Pack( 1, args );
	PyObject* copy = PyObject_CallObject( (PyObject*) &PositionWeightMatrixType, tuple );
	Py_DECREF( args );
	Py_DECREF( tuple );
	return copy;
	
	//return (PyObject*) self;
}

static PyObject* PositionWeightMatrix_reduce( PositionWeightMatrix* self )
{
	//returns tuple of ( type, constructor args )
	PyObject* args = PositionWeightMatrix_toList( self );
	PyObject* tuple = PyTuple_Pack( 2, &PositionWeightMatrixType, PyTuple_Pack( 1, args ) );
	//PyObject* copy = PyObject_CallObject( (PyObject*) &PositionWeightMatrixType, tuple );
	//Py_DECREF( args );
	//Py_DECREF( tuple );
	return tuple;
	
}

static PyObject* PositionWeightMatrix_nullScoreDistribution( PositionWeightMatrix* self, PyObject* args )
{
	PyObject* seq;
	int reps;
		
	if( !PyArg_ParseTuple( args, "Oi", &seq, &reps ) )
		return NULL;
	
	seq = PySequence_Fast( seq, "argument must be iterable" );
	if( !seq ) return NULL;
	
	PyObject** elems = PySequence_Fast_ITEMS( seq );
	size_t len = PySequence_Fast_GET_SIZE( seq );
	
	//allocate array to put elems in
	long* array = malloc( len * sizeof(long) );
	size_t i;
	for( i = 0 ; i < len ; ++i )
	{
		array[i] = PyInt_AsLong( elems[i] );
	}
	//check for error
	if( PyErr_Occurred() )
	{
		free( array );
		return NULL;
	}

	double** matrix = self->matrix;
	int n = self->n;
	int m = self->m;
	
	//allocate scores array
	size_t nmers = ( len - n + 1 );
	size_t lenScores = nmers*reps;
	double* scores = malloc( lenScores * sizeof( double ) ); 
	
	//for each replicate, shuffle and score the sequence
	int rep = 0;
	while( rep < reps )
	{
		shuffle( array, len, sizeof( long ) );
		scoreAllLong( matrix, n, m, array, len, scores, nmers * rep );
		++rep;
	}

	//sort the scores
	quicksort( scores, lenScores, sizeof( double ), &compareDoubleDescending );
	
	//free the long array
	free( array );
	
	//build a python list to return
	PyObject* list = PyList_New( lenScores );
	for( i = 0 ; i < lenScores ; ++i )
	{
		PyList_SetItem( list, i, PyFloat_FromDouble( scores[i] ) );
	}
	
	//free the scores array
	free( scores );
	
	return list;

}

static PyMemberDef PositionWeightMatrix_members[] = {
	{NULL} /* Sentinel */
};

static PyMethodDef PositionWeightMatrix_methods[] = {
	{ "score", (PyCFunction) PositionWeightMatrix_score, METH_VARARGS,
		"Scores the given sequence"
	},
	{ "scoreAll", (PyCFunction) PositionWeightMatrix_scoreAll, METH_VARARGS,
		"Scores all k-mers in the given sequence"
	},
	{"length", (PyCFunction) PositionWeightMatrix_length, METH_NOARGS,
		"Returns the length of this PWM"
	},
	{"nullScoreDistribution", (PyCFunction) PositionWeightMatrix_nullScoreDistribution, METH_VARARGS,
		"Generates a null score distribution by shuffling and scoring the sequence the specified number of times"
	},
	{"toList", (PyCFunction) PositionWeightMatrix_toList, METH_NOARGS,
		"Converts the PWM to a nested python list"
	},
	{"__copy__", (PyCFunction) PositionWeightMatrix_copy, METH_NOARGS,
		"Copies this PWM"
	},
	{"__deepcopy__", (PyCFunction) PositionWeightMatrix_deepcopy, METH_KEYWORDS,
		"Deeply copies this PWM"
	},
	{"__reduce__", (PyCFunction) PositionWeightMatrix_reduce, METH_NOARGS,
		"Provides pickling info for this class"
	},
	{NULL} /* Sentinel */
};

static PySequenceMethods PositionWeightMatrix_sequence_methods = {
	(lenfunc) PositionWeightMatrix_len, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	(ssizeargfunc) PositionWeightMatrix_getItem /*sq_item*/
};

static PyTypeObject PositionWeightMatrixType = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"PositionWeightMatrix.PositionWeightMatrix",             /*tp_name*/
    	sizeof(PositionWeightMatrix),             /*tp_basicsize*/
    	0,                         /*tp_itemsize*/
    	(destructor)PositionWeightMatrix_dealloc, /*tp_dealloc*/
    	0,                         /*tp_print*/
    	0,                         /*tp_getattr*/
    	0,                         /*tp_setattr*/
    	0,                         /*tp_compare*/
    	0,                         /*tp_repr*/
    	0,                         /*tp_as_number*/
    	&PositionWeightMatrix_sequence_methods,                         /*tp_as_sequence*/
    	0,                         /*tp_as_mapping*/
    	0,                         /*tp_hash */
    	0,                         /*tp_call*/
    	0,                         /*tp_str*/
    	0,                         /*tp_getattro*/
    	0,                         /*tp_setattro*/
    	0,                         /*tp_as_buffer*/
    	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    	"Position weight matrix",           /* tp_doc */
    	0,		               /* tp_traverse */
    	0,		               /* tp_clear */
    	0,		               /* tp_richcompare */
    	0,		               /* tp_weaklistoffset */
    	0,		               /* tp_iter */
    	0,		               /* tp_iternext */
    	PositionWeightMatrix_methods,             /* tp_methods */
    	PositionWeightMatrix_members,             /* tp_members */
    	0,                         /* tp_getset */
    	0,                         /* tp_base */
    	0,                         /* tp_dict */
    	0,                         /* tp_descr_get */
    	0,                         /* tp_descr_set */
    	0,                         /* tp_dictoffset */
    	(initproc)PositionWeightMatrix_init,      /* tp_init */
    	0,                         /* tp_alloc */
    	PositionWeightMatrix_new,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initPositionWeightMatrix(void) 
{
	PyObject* m;

    	if (PyType_Ready(&PositionWeightMatrixType) < 0)
    		return;

   	m = Py_InitModule3("PositionWeightMatrix", module_methods,
                       "Position weight matrix module.");

   	if (m == NULL)
   		return;

   	Py_INCREF(&PositionWeightMatrixType);
  	PyModule_AddObject(m, "PositionWeightMatrix", (PyObject *)&PositionWeightMatrixType);
}

