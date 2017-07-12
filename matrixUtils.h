#ifndef MATRIXUTILS_
#define MATRIXUTILS_

/*\brief Allocate space for a vector of doubles of length l
 * 
 * \param length Length of vector */
double *dvector(int length);

/*\brief Allocate space for a vector of ints of length l
 * 
 * \param length Length of vector */
int *ivector(int length);


/**\brief frees a vector created using dvector
 * 
 * \param v Vector to be freed */
void free_dvector(double *v);

/**\brief frees a vector created using ivector
 * 
 * \param v Vector to be freed */
void free_ivector(int *v);

/**\brief allocate integer matrix of size rows x cols
 * 
 * matrix is accessed via pointer to pointer, 
 * standard indexing is valid e.g. matrix[i][j]
 * 
 * Example usage:
 * int **matrix;
 * matrix=imatrix(3,3);
 * 
 * \param rows number of rows
 * \param cols number of columns */
int **imatrix(int rows, int cols);

/**\brief Frees a matrix created using imatrix
 * 
 * \param **a Matrix to be deallocated
 * \param rows Row size of a
 * \param cols Column size of a */
void free_imatrix(int **a, int rows, int cols);

/**
 * sum = vector + constant;
 * @param vector
 * @param constant
 * @param size
 * @param sum
 */
void add_int2ivec(int *vector, int constant, int size, int *sum);

void zeroIMatrix(int **matrix, int rows, int cols);

void printIVector(int *ivector, int n, char *msg);

void printDVector(double *dvector, int n, char *msg);
 
void printStr(char *msg);

void printIMatrix(int **imatrix, int nRow, int nCols);

double meanIVector(int arrSize, int *a);

int sumIVector(int arrSize, int *a);
/**
 * simple bubble sort algorithm for array of integers: sort in ascending order.
 *
 * @param a integer array
 * @param n array size
 * @param aIndices the indeces of the reordered array a
 */
void bubbleSrtIVect(int *a, int n, int *aIndices);

/**
 *  Create an array of doubles
 * @param array
 * @param arraySize
 */
void createRandDoubleArray(double *array, int arraySize);

/**
 * Initialize double array with 1.0, ..., arraySize
 * @param array
 * @param arraySize
 */
void createUniSpaceDoubleArray(double* array, int arraySize);

/**
 * Return a permutation matrix given the vector of permutations 
 * @param permuteVector
 * @param arraySize
 * @return 
 */
int **createPermutationMatrix(int *permuteVector, int arraySize);

/**@brief multiplies square matrices :
 *  [matrix3] = [matrix1] x [matrix2] 
 * 
 * @param **matrix1 Left input matrix
 * @param **matrix2 Right input matrix
 * @param **matrix3 Output matrix 
 * @param **n Matrix size
 */
void multiplySquareMatrix(double **matrix1, double **matrix2, double **matrix3, int n); 

/**
 * Reshape an integer matrix of nRows x nCols into a vector of length nRows*nCols
 * @param matrix
 * @param vector
 * @param nRows
 * @param nCols
 */
void reshapeIMatrix2vect(int **matrix, int *vector, int nRows, int nCols);

#endif
