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

#define MAXARRSIZE  40
#define LOAD_TRANSFER_CUTOFF 2

/**
 * Initialize each rank with an array of random size containing random doubles
 * @param xArr
 * @param resultsArr
 * @param arrSize
 */
void initRanks0(double **xArr, double **resultsArr, int *arrSize)
{

	*arrSize = rand() % MAXARRSIZE;
	*xArr = dvector(*arrSize);
	//createRandDoubleArray(xArr, *arrSize);
	createUniSpaceDoubleArray(*xArr, *arrSize);
	*resultsArr = dvector(*arrSize);
}

void initRanks(double **xArr, double **resultsArr, int *arrSize)
{

	*arrSize = rand() % MAXARRSIZE;
	*xArr = dvector(*arrSize);
	//createRandDoubleArray(xArr, *arrSize);
	createUniSpaceDoubleArray(*xArr, *arrSize);
	*resultsArr = dvector(*arrSize);
}

/**
 * Calculate the number of elements that should be transferred between ranks for equal
 * load between ranks.
 * To calculate the number of elements first sort the array loadSizes, then start
 * redistributing from both ends of the sorted array: i.e. start by sending some load from
 * the rank with the highest load to the rank with the lowest load. Then move inward 
 * until the loads are balanced. 
 *  
 * @param numProcs
 * @param initArraySizesAll
 * @param transferMatrix : t_ij is number of elements that rank i should receive from rank j
 * @param loadTransferCutoff : Only transfer load if above this value
 */
void calcRebalencingTransfers(int numProcs, int *initArraySizesAll, int *transferVectorAll)
{
	int loadAve; 	
	loadAve = (int) meanIVector(numProcs, initArraySizesAll);
	
	/*
	 * Sort the array to find the ranks with lowest to highest loads
	 */

	int *sortedIndices; // mapping from initArraySizeAll to the sorted array
	sortedIndices = ivector(numProcs);
	bubbleSrtIVect(initArraySizesAll, numProcs, sortedIndices);

	/*
	 * Find the difference between the initial load on each rank and the average load
	 */
	int *loadOffset; // difference between sorted array load and loadAve
	loadOffset = ivector(numProcs);
	add_int2ivec(initArraySizesAll, -loadAve, numProcs, loadOffset);
	
	/*
	 * Find how much load needs to be redistributed:
	 * Start from transferring from the highest loaded ranks to the lowest loaded ranks 
	 * and move inward.  
	 */
	int iDonor;
	int iRecep;
	int load2Trans;
	int **transferMatrixSorted;
	transferMatrixSorted = imatrix(numProcs, numProcs);
	zeroIMatrix(transferMatrixSorted, numProcs, numProcs);
	iDonor = numProcs -1;
	iRecep = 0;

	while (((initArraySizesAll[iRecep] - loadAve) < -LOAD_TRANSFER_CUTOFF) &&
		(iDonor > iRecep)) {
		while ((loadOffset[iRecep] < -LOAD_TRANSFER_CUTOFF) && (iDonor > iRecep)) {
			
			if (-loadOffset[iRecep] > loadOffset[iDonor]) {
				load2Trans = loadOffset[iDonor];
			} else {
				load2Trans = -loadOffset[iRecep];
			}
			transferMatrixSorted[iRecep][iDonor] = load2Trans;
			transferMatrixSorted[iDonor][iRecep] = -load2Trans;
			loadOffset[iRecep] += load2Trans;
			loadOffset[iDonor] -= load2Trans;
			
			if (-loadOffset[iRecep] > loadOffset[iDonor]) {
				iDonor--;
			}
		}
		iRecep++;
	}

	/*
	 * Reorder transferMatrix due to initial sorting
	 */	
	int **transferMatrix;
	
	transferMatrix = imatrix(numProcs, numProcs);
	zeroIMatrix(transferMatrix, numProcs, numProcs);
	int iRow, iRow_actual;
	int iCol, iCol_actual;
	for (iRow = 0; iRow < numProcs; iRow++) {
		for (iCol = 0; iCol < numProcs; iCol++) {
			iRow_actual = sortedIndices[iRow];
			iCol_actual = sortedIndices[iCol];
			transferMatrix[iRow_actual][iCol_actual] = transferMatrixSorted[iRow][iCol];
		}
	}
	
	printf("transferMatrix: \n"); printIMatrix(transferMatrix, numProcs, numProcs);
	reshapeIMatrix2vect(transferMatrix, transferVectorAll, numProcs, numProcs);
	//printf("transferVectorAll:  "); 
/*
	char vectMsg[40];
	strcpy(vectMsg, "transferVectorAll:   " );
	printIVector(transferVectorAll, numProcs * numProcs, vectMsg);
*/
	free_imatrix(transferMatrix, numProcs, numProcs);
	
	free_imatrix(transferMatrixSorted, numProcs, numProcs);
	free_ivector(loadOffset);
	free_ivector(sortedIndices);
}

enum rank_type {
	RECEIPIENT = -1,
	NOCOMM = 0,
	DONOR = 1
};

struct transfer_info {
	int initArraySize;
	int totalLoad2transfer;
	int nTransfers;
	enum rank_type rankType;
	int *ranks2comm;
	int *arraySizes2comm;
};

void printRankTransferInfo(struct transfer_info transInfo, int myRank)
{
	char prntMsg[1024];
	char tmpStr[1000];
	int iRank;
	
	switch (transInfo.rankType) {
	case RECEIPIENT:

		sprintf(tmpStr, "Rank %d receives", myRank);
		strcpy(prntMsg, tmpStr);

		for (iRank = 0; iRank < transInfo.nTransfers; iRank++) {
			sprintf(tmpStr, " (%d->) %d ", transInfo.ranks2comm[iRank]
				, transInfo.arraySizes2comm[iRank]);
			strcat(prntMsg, tmpStr);
		}
		sprintf(tmpStr, " : total = %d", transInfo.totalLoad2transfer);
		strcat(prntMsg, tmpStr);
		printf("%s \n", prntMsg);

		break;
	case DONOR:
		sprintf(tmpStr, "Rank %d sends", myRank);
		strcpy(prntMsg, tmpStr);
		for (iRank = 0; iRank < transInfo.nTransfers; iRank++) {
			sprintf(tmpStr, " %d (->%d)", -transInfo.arraySizes2comm[iRank],
				transInfo.ranks2comm[iRank]);
			strcat(prntMsg, tmpStr);
		}
		sprintf(tmpStr, " : total = %d", -transInfo.totalLoad2transfer);		
		strcat(prntMsg, tmpStr);
		printf("%s \n", prntMsg);
		break;
	case NOCOMM:
		printf("Rank %d has just about the right amount of work to do!\n", myRank);
	}

}

void getRankTransferInfo(int *transferArr, int numProcs, struct transfer_info *transInfo)
{
	
	transInfo->nTransfers = 0;
	transInfo->totalLoad2transfer = sumIVector(numProcs, transferArr);
	
	/*
	 * What type of rank is this?
	 */
	if (transInfo->totalLoad2transfer < 0) {
		transInfo->rankType = DONOR;
	} else if (transInfo->totalLoad2transfer == 0) {
		transInfo->rankType = NOCOMM;
	} else {
		transInfo->rankType = RECEIPIENT;
	}

	int iRank;
	for (iRank = 0; iRank < numProcs; iRank++) {
		if (transferArr[iRank] !=0){
			transInfo->nTransfers ++;
		}
	}
	
	transInfo->ranks2comm = ivector(transInfo->nTransfers);
	transInfo->arraySizes2comm = ivector(transInfo->nTransfers);
	int count = 0;
	for (iRank = 0; iRank < numProcs; iRank++) {
		if (transferArr[iRank] !=0){
			transInfo->ranks2comm[count] = iRank;
			transInfo->arraySizes2comm[count] = transferArr[iRank];
			count ++;
		}
	}
}

void calc_local(double *x, double *sin_x, int nArray){
	int i;
	for (i = 0; i < nArray; i++) {
		sin_x[i] = sin(x[i]);
	}
}

void trasnfer_data2calc(struct transfer_info myTransferInfo, double **extraLoad, 
	double *localArr, MPI_Request *data2calcReqs, int nArr0myRank, int myRank)
{

	int iRank;
	int bufferIndex;
	
	if (myTransferInfo.rankType == RECEIPIENT) {
		*extraLoad = dvector(myTransferInfo.totalLoad2transfer);
	}
	
	switch (myTransferInfo.rankType) {
	case RECEIPIENT:
	{
		for (iRank = 0; iRank < myTransferInfo.nTransfers; iRank++) {

			bufferIndex = sumIVector(iRank, myTransferInfo.arraySizes2comm);
/*
			printf("R %d: from %d : extraLoad[bufferIndex] %f",
				myRank,myTransferInfo.ranks2comm[iRank]);
*/
			MPI_Irecv(&((*extraLoad)[bufferIndex]),
				myTransferInfo.arraySizes2comm[iRank],
				MPI_DOUBLE,
				myTransferInfo.ranks2comm[iRank],
				0,
				MPI_COMM_WORLD,
				&data2calcReqs[iRank]);

		}
		break;
	}
	case DONOR:
	{
		double *load2send;
		int overloadStartIndex;
		overloadStartIndex = nArr0myRank + myTransferInfo.totalLoad2transfer;
		load2send = &localArr[overloadStartIndex];

		for (iRank = 0; iRank < myTransferInfo.nTransfers; iRank++) {
			bufferIndex = -sumIVector(iRank, myTransferInfo.arraySizes2comm);
			
			char prntMsg[100];
			sprintf(prntMsg,"R %d:load2send -> %d, (bufferIndex : %d )",
				myRank, 
				myTransferInfo.ranks2comm[iRank], bufferIndex );
			
			printDVector(&load2send[bufferIndex],-myTransferInfo.arraySizes2comm[iRank],prntMsg);
			
			MPI_Isend(&load2send[bufferIndex],
				-myTransferInfo.arraySizes2comm[iRank],
				MPI_DOUBLE,
				myTransferInfo.ranks2comm[iRank],
				0,
				MPI_COMM_WORLD,
				&data2calcReqs[iRank]);

		}
		break;
	}
	case NOCOMM:
		break;
	}
}

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
	char msgIvector[40];
	strcpy(msgIvector, " test print: "); 
	
	printIVector(vect,2*3, msgIvector);
	printStr(msgIvector);
	printf("*****************************\n");

	free_ivector(vect);
	
}



/*
int main(int argc, char**argv)
{
	test_reshapeIMatrix();
	exit(0);
}
*/
int main(int argc, char**argv){
	int myRank;
	int numProcs;

	/* initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	/*
	 * Initialize loads on each rank
	 */

	int nArr0myRank; // initial size of the array
	double *localDataArr, *sinAnglesArr;
	srand(time(NULL) + myRank);  // To get different random numbers on each processor
	//int	givenArrsizes[5] = {31,   32,    3,   27,   24};
	initRanks(&localDataArr, &sinAnglesArr, &nArr0myRank);
	//nArr0myRank = givenArrsizes[myRank];
	//initRanks(&localDataArr, &sinAnglesArr, nArr0myRank);

	//printf("Hello from proc %d! My initial load is %d \n", myRank, nArr0myRank);
	//MPI_Barrier(MPI_COMM_WORLD);
	
	/*
	 * Gather initial load sizes from all ranks 
	 */
	int *initArraySizesAll = NULL;
	if (myRank == 0) 
		initArraySizesAll = ivector(numProcs);
	
	MPI_Gather(&nArr0myRank, 1, MPI_INT, initArraySizesAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (myRank == 0) {
		char msgIvector[40];
		strcpy(msgIvector, " initial loads"); 
		printIVector(initArraySizesAll, numProcs, msgIvector);
	}
	
	/*
	 * Calculate the number of transfers required between ranks for equal balancing
	 */
	
	int *transferVectorAll = NULL; 
	if (myRank == 0) {
		transferVectorAll = ivector(numProcs * numProcs);		
		calcRebalencingTransfers(numProcs, initArraySizesAll, transferVectorAll);
	}

	/*
	 * Transfer each row of the transferMatrix to a processor
	 */
	int *myTrasferArr = ivector(numProcs);
	MPI_Scatter(transferVectorAll, numProcs , MPI_INT,
		myTrasferArr,numProcs, MPI_INT, 0, MPI_COMM_WORLD);
	free_ivector(transferVectorAll);		
	/*
	 * Redistribute arrays before starting to work
	 */
	
	struct transfer_info myTransferInfo; 
	getRankTransferInfo(myTrasferArr, numProcs, &myTransferInfo);
	printRankTransferInfo(myTransferInfo, myRank);
	
	/*
	 * Loop through the transfer array and send/receive data if necessary.
	 */
	

	if (myTransferInfo.rankType != NOCOMM) {
		MPI_Request *data2calcReqs;
		MPI_Status *data2calcStats;
		double *extraLoad;
		
		data2calcReqs = (MPI_Request *) malloc(myTransferInfo.nTransfers * sizeof(MPI_Request));
		data2calcStats = (MPI_Status *) malloc(myTransferInfo.nTransfers * sizeof(MPI_Status));

		if (myTransferInfo.rankType == RECEIPIENT) {
			char vectMssg[100];
			sprintf(vectMssg, "ExtraLoad buffer before receive %d (size(%d)) :", myRank,
				myTransferInfo.totalLoad2transfer);
			printDVector(extraLoad, myTransferInfo.totalLoad2transfer, vectMssg);
		}

		trasnfer_data2calc(myTransferInfo, &extraLoad, localDataArr, data2calcReqs,
			nArr0myRank, myRank);


		MPI_Waitall(myTransferInfo.nTransfers, data2calcReqs, data2calcStats);

		if (myTransferInfo.rankType == RECEIPIENT) {
			char vectMssg[200];
			sprintf(vectMssg,"transfered load to %d : ", myRank);
			printDVector(extraLoad, myTransferInfo.totalLoad2transfer, vectMssg);
		}

		free(data2calcReqs);
		free(data2calcStats);
		if (myTransferInfo.rankType == RECEIPIENT) {
			free_dvector(extraLoad);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * Clean up
	 */

	free_ivector(myTrasferArr);
	free_dvector(localDataArr);
	free_dvector(sinAnglesArr);
	free_ivector(myTransferInfo.arraySizes2comm);
	free_ivector(myTransferInfo.ranks2comm);


	if (myRank == 0) {
		free_ivector(initArraySizesAll);
	}

	MPI_Finalize(); /* cleanup MPI */
	exit(0);
}


