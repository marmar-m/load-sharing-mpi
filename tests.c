/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: Marmar
 *
 * Created on July 10, 2017, 8:54 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "matrixUtils.h"

void test_reshapeIMatrix();

void test_reshapeIMatrix()
{
	int **testMat;
	int *vect;
	testMat = imatrix(2,3);
	int i,j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 3; j++) {
			testMat[i][j] = i*3+j+1;
		}
	}
	printf("*****************************\n");
	printIMatrix(testMat,2,3);
	printf("*****************************\n");
	vect = ivector(6);
	reshapeIMatrix2vect(testMat,vect,2,3);
	
	printf("*****************************\n");
	printIVector(vect,2*3);
	printf("*****************************\n");

	free_ivector(vect);
	
}


