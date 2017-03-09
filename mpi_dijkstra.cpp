/*
 * This is code skeleton for COMP5112 assignment1
 * Compile: mpic++ -o mpi_dijkstra mpi_dijkstra.cpp
 * Run: mpiexec -n <number of processes> mpi_dijkstra <input file>, you will find the output in 'output.txt' file
 */


#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <cstring>
#include <algorithm>
#include "mpi.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;


/*
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and one matrix dimension convert(2D->1D) function
 */
namespace utils {
    int N; //number of vertices
    int *mat; // the adjacency matrix

    /*
     * convert 2-dimension coordinate to 1-dimension
     */
    int convert_dimension_2D_1D(int x, int y) {
        return x * N + y;
    }

    int read_file(string filename) {
        std::ifstream inputf(filename, std::ifstream::in);
        inputf >> N;
        assert(N < (1024 * 1024 *
                    20)); // input matrix should be smaller than 20MB * 20MB (400MB, we don't have two much memory for multi-processors)
        mat = (int *) malloc(N * N * sizeof(int));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                inputf >> mat[convert_dimension_2D_1D(i, j)];
            }

        return 0;
    }

    string format_path(int i, int *pred) {
        string out("");
        int current_vertex = i;
        while (current_vertex != 0) {
            string s = std::to_string(current_vertex);
            std::reverse(s.begin(), s.end());
            out = out + s + ">-";
            current_vertex = pred[current_vertex];
        }
        out = out + std::to_string(0);
        std::reverse(out.begin(), out.end());
        return out;
    }

    int print_result(int *dist, int *pred) {
        std::ofstream outputf("output.txt", std::ofstream::out);
        outputf << dist[0];
        for (int i = 1; i < N; i++) {
            outputf << " " << dist[i];
        }
        for (int i = 0; i < N; i++) {
            outputf << "\n";
            if (dist[i] >= 1000000) {
                outputf << "NO PATH";
            } else {
                outputf << format_path(i, pred);
            }
        }
        outputf << endl;
        return 0;
    }
}//namespace utils
// you may add some helper functions here.

//get min Distance from distance list
int getMin(int *dist, bool *visit, int loc_n, int loc_N, int p, int my_rank)
{
    int minIndex = -1;
    //find the first one that is unvisited, if all visited, return -1
    for(int i = 0; i < loc_n; i++)
    {
        if (!visit[i])
        {
            minIndex = i;
            break;
        }
    }
    if (minIndex == -1)
    {
        return -1;
    }
    for(int i = 0; i < loc_n; i++)
    {
        if (!visit[i] && dist[i] < dist[minIndex])
        {
            minIndex = i;
        }
    }
    return minIndex;
}

void dijkstra(int my_rank, int N, int p, MPI_Comm comm, int *mat, int *all_dist, int *all_pred) {

    //------your code starts from here------
    int loc_N; // I need a local copy for N
    int loc_n; //how many vertices I need to process.
    int *loc_mat; //local matrix
    int *loc_dist; //local distance
    int *loc_pred; //local predecessor
    bool *loc_visit; //local visit record array

    //step 1: broadcast N
    loc_N = N;
    //step 2: find loc_n
    if (my_rank == p-1)
    {
        loc_n = N - N/p * p + N/p; 
    } 
    else 
    {
        loc_n = N/p;
    }
    //step 3: allocate local memory
    loc_mat = (int *) malloc(loc_N*loc_n*sizeof(int));
    loc_dist = (int *) malloc(loc_n*sizeof(int));
    loc_pred = (int *) malloc(loc_n*sizeof(int));
    loc_visit = (bool *) malloc(loc_n*sizeof(bool));
    //step 4: broadcast matrix mat
    for (int i = 0; i < loc_N; i++)
    {
        for (int j =0; j < loc_n; j++)
        {
            loc_mat[i * loc_n + j] = mat[i * loc_N + j + my_rank* (N/p)];
        }
    }
    for (int i = 0; i < loc_n; i++)
    {
        loc_dist[i] = 1000000;
        loc_visit[i] = 0;
    }
    if (my_rank == 0)
    {
        loc_dist[0] = 0;
    }

    //step 4: dijkstra algorithm

    int Infinity = 1000001;


for(int iter = 0; iter < 5000; iter++)
{   
    int minGlobal = 0;
    int minDistGlobal = 0;
    int minRank = 0;
    int loc_minIndex = getMin(loc_dist, loc_visit, loc_n, loc_N, p, my_rank);

    //Here we get the local Min Distance and It's Index...Send it out....
    if (my_rank != 0)
    {
        MPI_Send(&loc_minIndex, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        if (loc_minIndex >= 0)
        {
            MPI_Send(&loc_dist[loc_minIndex], 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        } else {
            MPI_Send(&Infinity, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
    } else {
        int *minIndex;
        minIndex = (int *)malloc(p*sizeof(int));
        minIndex[0] = loc_minIndex;
        int *minDist;
        minDist = (int *)malloc(p*sizeof(int));
        if (loc_minIndex >= 0)
        {
            minDist[0] = loc_dist[loc_minIndex];
        } else {
            minDist[0] = Infinity;
        }
        for (int i = 1; i < p; i++)
        {
            MPI_Recv(&minIndex[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&minDist[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int i = 0; i < p; i++)
        {
            if(minDist[i] < minDist[minGlobal])
            {
                minGlobal = i;
                minRank = i;
            }
        }
        minDistGlobal = minDist[minGlobal];
        minGlobal = minGlobal * (loc_N/p) + minIndex[minGlobal];

        for (int i = 1; i < p; i++)
        {
            MPI_Send(&minGlobal, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
            MPI_Send(&minDistGlobal, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            MPI_Send(&minRank, 1, MPI_INT, i, 4, MPI_COMM_WORLD);
        }
    }

    if (my_rank != 0) {
        MPI_Recv(&minGlobal, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&minDistGlobal, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&minRank, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (int i = 0; i < loc_n; i++)
    {
        if(!loc_visit[i])
        {
        if (minDistGlobal + loc_mat[minGlobal * loc_n + i] < loc_dist[i])
        {
            loc_dist[i] = minDistGlobal + loc_mat[minGlobal * loc_n + i];
            loc_pred[i] = minGlobal;
        }
    }
    }

    if (my_rank == minRank)
    {
        int n = minGlobal - minRank * (loc_N/p);
        loc_visit[n] = 1;
    }
}
       //step 5: retrieve results back
    //Hint: use MPI_Gather(or MPI_Gatherv) function
    int* rcounts;
    int* displs;
    rcounts = (int *) malloc(p*sizeof(int));
    displs = (int *) malloc(p*sizeof(int));
 
    int temp = loc_N/p;
    for (int i = 0; i < p; i++)
    {
        rcounts[i] = temp;
        displs[i] = i * temp;
    }
    rcounts[p-1] = N - N/p * p + N/p;

    //MPI_Gather(loc_dist, loc_n, MPI_INT, all_dist, loc_n, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Gather(loc_pred, loc_n, MPI_INT, all_pred, loc_n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(loc_dist, loc_n, MPI_INT, all_dist, rcounts, displs, MPI_INT, 0 , MPI_COMM_WORLD);
    MPI_Gatherv(loc_pred, loc_n, MPI_INT, all_pred, rcounts, displs, MPI_INT, 0 , MPI_COMM_WORLD);
    //step 6: remember to free memory
    free(loc_dist);
    free(loc_mat);
    free(loc_visit);

    //------end of your code------
}


int main(int argc, char **argv) {
    assert(argc > 1 && "Input file was not found!");
    string filename = argv[1];
    assert(utils::read_file(filename) == 0);

    //`all_dist` stores the distances and `all_pred` stores the predecessors
    int *all_dist;
    int *all_pred;
    all_dist = (int *) calloc(utils::N, sizeof(int));
    all_pred = (int *) calloc(utils::N, sizeof(int));

    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int p;//number of processors
    int my_rank;//my global rank

    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    dijkstra(my_rank, utils::N, p, comm, utils::mat, all_dist, all_pred);

    if (my_rank == 0)
        utils::print_result(all_dist, all_pred);
    MPI_Finalize();


    free(utils::mat);
    free(all_dist);
    free(all_pred);

    return 0;
}
