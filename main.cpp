#include <cassert>
#include <cstdint>
#include <cstdio>

#include <cmath>

#include <algorithm>
#include <utility>

#include <mpi.h>

#include "std_LCG_PLE63.hpp"
#include "uniform_distribution.hpp"

#define DARTS 500000     /* number of throws at dartboard */
#define ROUNDS 100      /* number of times "darts" is iterated */
#define MASTER 0        /* task ID of master task */

using ugen_t = std::linear_congruential_engine<uint64_t, 2806196910506780709ULL, 1ULL, (1ULL<<63ULL)>;

double dboard(int darts, ugen_t& ugen);

int main (int argc, char *argv[])
{
    int	taskid,	        /* task ID - also used as seed number */
	numtasks;       /* number of tasks */
    MPI_Status status;

    /* Obtain number of tasks and task ID */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    printf ("MPI task %d has started...\n", taskid);
    
    ugen_t ugen;
    ugen.discard(taskid*1000000000);
    
    double pi = 0.0;
    double avepi = 0.0;
    for (int i = 0; i != ROUNDS; ++i)
    {
        /* All tasks calculate pi using dartboard algorithm */
        double homepi = dboard(DARTS, ugen);

        /* Workers send homepi to master */
        /* - Message type will be set to the iteration count */
        if (taskid != MASTER)
        {
            int mtype = i;
            int rc = MPI_Send(&homepi, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
            if (rc != MPI_SUCCESS)
                printf("%d: Send failure on round %d\n", taskid, mtype);
        } 
        else
        {
            int mtype = i;
            double pisum = 0.0;
            for (int n = 1; n < numtasks; n++)
            {
                double pirecv;
                int rc = MPI_Recv(&pirecv, 1, MPI_DOUBLE, MPI_ANY_SOURCE, mtype, MPI_COMM_WORLD, &status);
                if (rc != MPI_SUCCESS) 
                    printf("%d: Receive failure on round %d\n", taskid, mtype);
                /* keep running total of pi */
                pisum = pisum + pirecv;
            }
            /* Master calculates the average value of pi for this iteration */
            pi = (pisum + homepi)/numtasks;
            avepi = ((avepi * i) + pi)/(i + 1); 
            printf("   After %8d throws, average value of pi = %10.8f\n", (DARTS * (i + 1)),avepi);
        }
    } 
    
    if (taskid == MASTER) 
        printf ("\nReal value of PI: 3.1415926535897 \n");

    MPI_Finalize();
    
    return 0;
}


static inline double squared(double x)
{
    return x*x;
}

double dboard(int darts, ugen_t& ugen)
{
    std::uniform_distribution<double> rng;
    int score = 0;

    double r;
    for (int n = 1; n <= darts; n++)
    {
        /* generate random numbers for x and y coordinates */
        r = rng(ugen);
        double x_coord = (2.0 * r) - 1.0;
        r = rng(ugen);
        double y_coord = (2.0 * r) - 1.0;

        /* if dart lands in circle, increment score */
        if ((squared(x_coord) + squared(y_coord)) <= 1.0)
           score++;
    }

    /* calculate pi */
    double pi = 4.0 * (double)score/(double)darts;
    
    return pi;
} 

