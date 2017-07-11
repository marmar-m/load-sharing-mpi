/*Utilities for matrices and vectors */
#include <stdio.h>
#include <stdlib.h>

/** Macro for swapping int variables */
#define SWAPI(a,b)   { int t; t = a; a = b; b = t; }  

double *dvector(int length)
{
	double *v;
	v=(double *) malloc(length*sizeof(double));
	if(v==NULL){
		printf("dvector allocation error of size %d\n", length);
		exit(1);
	}
	return v;
}

int *ivector(int length)
{
	int *v;
	v=(int*) malloc(length*sizeof(int));
	if(v==NULL){
		printf("ivector allocation error of size %d\n", length);
		exit(1);
	}
	return v;
}

void free_dvector(double *v)
{
	free(v);
}

void free_ivector(int *v)
{
	free(v);
}

int **imatrix(int rows, int cols)
{
	/* Use contiguous memory--row major (i*cols + j) */
	int i;
	int **m;
	int *v;

	/* allocate memory in bulk */
	v = (int *) malloc(rows*cols*sizeof(int));
	if(v == NULL){
		printf("imatrix allocation error (%d, %d)\n", rows, cols);
		exit(1);
	}

	/* allocate pointers */
	m=(int **) malloc(rows*sizeof(int *));
	if(m == NULL){
			printf("imatrix allocation error (%d, %d)\n", rows, cols);
			exit(1);
	}
	/* now point m to appropriate places in v */
	for (i = 0; i < rows; i++) {
		m[i] = &v[i * cols];
	}

	return m;
}

int **imatrix0(int n, int m)
{	
	int *data = (int *)malloc(n*m*sizeof(int));
    int **array = (int **)malloc(n*sizeof(int *));
    int i;
	for ( i=0; i<n; i++){
        array[i] = &(data[i*m]);
	}
    return array;
}

void free_imatrix(int **a, int rows, int cols)
{
	free(*a);
	free(a);
}


void add_int2ivec(int *vector, int constant, int size, int *sum)
{
    int i;
    for (i = 0; i < size; i++) {
        sum[i] = vector[i] + constant;
    }
}

void zeroIMatrix(int **matrix, int rows, int cols) {
    int row, column;
    for (row = 0; row < rows; row++) {
        for (column = 0; column < cols; column++) {
            matrix[row][column] = 0.0;
        }
    }
}

void printIVector(int *ivector, int n)
{
	int i;
	for (i = 0; i < n; i++) {
		printf("%5d ", ivector[i]);
	}
	printf("\n");
}


void printDVector(double *dvector, int n)
{
	int i;
	for (i = 0; i < n; i++) {
		printf("%3.2f ", dvector[i]);
	}
	printf("\n");
}


void printIMatrix(int **imatrix, int nRows, int nCols)
{
	int row, column;
	for (row = 0; row < nRows; row++) {
		for (column = 0; column < nCols; column++) {
			printf("%5d ", imatrix[row][column]);
		}
		printf("\n");
	}
	printf("\n");
}

double meanIVector(int arrSize, int *a) {
    int sum=0, i;
    for(i=0; i<arrSize; i++)
        sum+=a[i];
    return((double )sum/arrSize);
}


int sumIVector(int arrSize, int *a) {
    int sum=0, i;
    for(i=0; i<arrSize; i++)
        sum+=a[i];
    return sum;
}

void bubbleSrtIVect(int *a, int n, int *aIndices) {
    int i, j;

    for (i = 0; i < n; i++) {
        aIndices[i] = i;
    }


    for (i = 0; i < n; i++) { // Make a pass through the array for each element
        for (j = 1; j < (n - i); j++) {// Go through the array beginning to end
            if (a[j - 1] > a[j]) { // If the the first number is greater, swap it 
                SWAPI(a[j - 1], a[j]);
                SWAPI(aIndices[j - 1], aIndices[j]);
            }
        }
    }
}

void createRandDoubleArray(double *array, int arraySize)
{
	int i;
	for (i = 0; i < arraySize; i++) {
		array[i] = rand();
	}
}

void createUniSpaceDoubleArray(double* array, int arraySize)
{
	int i;
	for (i = 0; i < arraySize; i++) {
		array[i] = (double)(i + 1);
	}

}

int **createPermutationMatrix(int *permuteVector, int arraySize){
	
	int **permuteMatrix;
	permuteMatrix = imatrix(arraySize, arraySize);
	zeroIMatrix(permuteMatrix, arraySize, arraySize);
	int i;
	for (i = 0; i < arraySize; i++) {
		permuteMatrix[permuteVector[i]][i] = 1;
	}
	return permuteMatrix;
}

void multiplySquareMatrix(double **matrix1, double **matrix2, double **matrix3, int n)
{
	int row, column, rowloop;
	for (row = 0; row < n; row++) {
		for (column = 0; column < n; column++) {
			matrix3[row][column] = 0.0;
			for (rowloop = 0; rowloop < n; rowloop++) {
				matrix3[row][column] += matrix1[row][rowloop] * matrix2[rowloop][column];
			}
		}
	}
}

void reshapeIMatrix2vect(int **matrix, int *vector, int nRows, int nCols){
	int iRow;
	int iCol;
	for (iRow = 0; iRow < nRows; iRow++) {
		for (iCol = 0; iCol < nCols; iCol++) {
			vector[iRow * nCols + iCol] = matrix[iRow][iCol];
		}
	}
}