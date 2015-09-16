#include <cassert>
#include <cstdint>
#include <cstdio>

#include <cmath>

#include <algorithm>
#include <utility>

#include <mpi.h>

#include "std_LCG_PLE63.hpp"
#include "uniform_distribution.hpp"

#define  MASTER		0

int main (int argc, char *argv[])
{
    int   numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);
    
    printf ("Hello from task %d on %s!\n", taskid, hostname);
    
    if (taskid == MASTER)
        printf("MASTER: Number of MPI tasks is: %d\n",numtasks);
        
    MPI_Finalize();

    return 0;
}

