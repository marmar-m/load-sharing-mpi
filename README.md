A Parallel Example with MPI
===========================

This is a C program with MPI to load balance a computation

1.  N MPI ranks
2.  Each MPI rank is initialized with a array of random size with random angles to compute ​sin(x)
3.  The program determines how to rebalance the load between the ranks and initializes non-blocking MPI send/receive commands to initiate the process of data transfer to rebalance load
4.  Call ​sin(x)​ on each rank with local data
5.  Wait for rebalance transfers to finish and perform calculations on newly received data, if any
6.  Send back the results to their original rank

Running the program
-------------------

-   The program requires `mpicc.mpich2` for compiling. After compiling with:

    ``` example
    > make
    ```

    An executable file (`loadsharingexample`) is generated under `dist/Debug/GNU-Linux/`.
-   The executable does not require inputs.

    ``` example
    mpirun.mpich -n 8 ./dist/Debug/GNU-Linux/loadsharingexample  
    ```

    The following sample results are generated

    ``` example
    initial loads: 1816 3150 4283 5375 9096 2571  304 8434  
    transferMatrix: 
    0     0     0     0   644     0     0  1918 
    0     0     0   897     0     0     0   331 
    0     0     0    95     0     0     0     0 
    0  -897   -95     0     0     0     0     0 
    -644     0     0     0     0     0 -4074     0 
    0     0     0     0     0     0     0  1807 
    0     0     0     0  4074     0     0     0 
    -1918  -331     0     0     0 -1807     0     0 

    Rank 0 receives (4->) 644  (7->) 1918  : total = 2562 
    Rank 0 : done with calculation of original array.
    Rank 0 : done with calculation of received data (2562) .
    Rank 0 : done with sending back data.
    Rank 1 receives (3->) 897  (7->) 331  : total = 1228 
    Rank 1 : done with calculation of original array.
    Rank 1 : done with calculation of received data (1228) .
    Rank 1 : done with sending back data.
    Rank 2 receives (3->) 95  : total = 95 
    Rank 2 : done with calculation of original array.
    Rank 2 : done with calculation of received data (95) .
    Rank 2 : done with sending back data.
    Rank 3 sends 897 (->1) 95 (->2) : total = 992 
    Rank 3 : done with calculation of original array.
    Rank 3 : done with sending back data.
    Rank 4 sends 644 (->0) 4074 (->6) : total = 4718 
    Rank 4 : done with calculation of original array.
    Rank 4 : done with sending back data.
    Rank 5 receives (7->) 1807  : total = 1807 
    Rank 5 : done with calculation of original array.
    Rank 5 : done with calculation of received data (1807) .
    Rank 5 : done with sending back data.
    Rank 6 receives (4->) 4074  : total = 4074 
    Rank 6 : done with calculation of original array.
    Rank 6 : done with calculation of received data (4074) .
    Rank 6 : done with sending back data.
    Rank 7 sends 1918 (->0) 331 (->1) 1807 (->5) : total = 4056 
    Rank 7 : done with calculation of original array.
    Rank 7 : done with sending back data.
    ```
