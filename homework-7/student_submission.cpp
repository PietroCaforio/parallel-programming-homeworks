//
// Created by Vincent Bode on 08/07/2020.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "life.h"
#include "Utility.h"
#include "VideoOutput.h"
#include <mpi.h>
#include <iostream>
#include <random>
#include <thread>
#include "Utility.h"

std::minstd_rand my_randomEngine;
uint_fast32_t my_cacheValue;
uint_fast32_t my_bitMask = 0;

void my_seedGenerator(unsigned long long seed)
{
    my_randomEngine = std::minstd_rand(seed);
}

inline bool generateBit()
{
    if (!my_bitMask)
    {
        my_cacheValue = my_randomEngine();
        my_bitMask = 1;
    }
    bool value = my_cacheValue & my_bitMask;
    my_bitMask = my_bitMask << 1;
    return value;
}

/*
  Apply the game of life rules on a Torus --> grid contains shadow rows and columns
  to simplify application of rules i.e. grid actually ranges from grid [ 1.. height - 2 ][ 1 .. width - 2]
*/
void evolve(ProblemData &problemData, int rank, int size)
{
    int workload_rows = GRID_SIZE / size;
    int workload_rest = GRID_SIZE % size;

    auto &grid = *problemData.readGrid;
    auto &writeGrid = *problemData.writeGrid;

    int flag = ((rank + 1) * workload_rows);

    if (rank == size - 1)
    {
        flag += workload_rest - 1;
    }

    // For each cell
    for (int i = rank == 0 ? 1 : rank * workload_rows; i < flag; i++)
    {
        for (int j = 1; j < GRID_SIZE - 1; j++)
        {
            // Calculate the number of neighbors
            int sum = grid[i - 1][j - 1] + grid[i - 1][j] + grid[i - 1][j + 1] +
                      grid[i][j - 1] + grid[i][j + 1] +
                      grid[i + 1][j - 1] + grid[i + 1][j] + grid[i + 1][j + 1];

            if (!grid[i][j])
            {
                // If a cell is dead, it can start living by reproduction or stay dead
                if (sum == 3)
                {
                    // reproduction
                    writeGrid[i][j] = true;
                }
                else
                {
                    writeGrid[i][j] = false;
                }
            }
            else
            {
                // If a cell is alive, it can stay alive or die through under/overpopulation
                if (sum == 2 || sum == 3)
                {
                    // stays alive
                    writeGrid[i][j] = true;
                }
                else
                {
                    // dies due to under or overpopulation
                    writeGrid[i][j] = false;
                }
            }
        }
    }
}

/*
  Copies data from the inner part of the grid to
  shadow (padding) rows and columns to transform the grid into a torus.
*/
void copy_edges(bool (&grid)[GRID_SIZE][GRID_SIZE], int rank, int size)
{
    int workload_rows = GRID_SIZE / size;
    int workload_rest = GRID_SIZE % size;

    int prev_proc = (rank - 1 + size) % size;
    int next_proc = (rank + 1) % size;

    // printf("PROC: %d  previous %d and next %d\n", rank, prev_proc, next_proc);

    // printf("Params: %d %d\n", workload_rows, workload_rest);

    if (rank != size - 1)
    {
        MPI_Request req[4];
        // Copy data to the boundaries
        for (int i = rank == 0 ? 1 : rank * workload_rows; i < ((rank + 1) * workload_rows); i++)
        {
            // join rows together
            grid[i][0] = grid[i][GRID_SIZE - 2];
            grid[i][GRID_SIZE - 1] = grid[i][1];
        }

        int david = rank == 0 ? 1 : 0;
        // send top row to previous process
        MPI_Isend(grid[rank * workload_rows + david], GRID_SIZE, MPI_CXX_BOOL, prev_proc, 0, MPI_COMM_WORLD, &req[0]);
        // printf("PROC %d, sending to %d\n", rank, prev_proc);

        // recv on bottom row from next
        // printf("PROC %d, recving to %d\n", rank, next_proc);
        MPI_Irecv(grid[(rank + 1) * workload_rows], GRID_SIZE, MPI_CXX_BOOL, next_proc, 0,
                  MPI_COMM_WORLD, &req[1]);

        // send bottom row to next process
        MPI_Isend(grid[(rank + 1) * workload_rows - 1], GRID_SIZE, MPI_CXX_BOOL, next_proc, 0, MPI_COMM_WORLD, &req[2]);
        // printf("PROC %d, sending to %d\n", rank, next_proc);

        // printf("PROC %d, recving to %d\n", rank, prev_proc);
        // recv on top row from previous process
        MPI_Irecv(grid[rank * workload_rows + david - 1], GRID_SIZE, MPI_CXX_BOOL, prev_proc, 0,
                  MPI_COMM_WORLD, &req[3]);

        MPI_Waitall(4, req, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Request req[4];

        // Copy data to the boundaries
        for (int i = rank * workload_rows; i < ((rank + 1) * workload_rows + workload_rest) - 1; i++)
        {
            // join rows together
            grid[i][0] = grid[i][GRID_SIZE - 2];
            grid[i][GRID_SIZE - 1] = grid[i][1];
        }

        // printf("PROC %d, recving to %d\n", rank, next_proc);
        // recv bottom padding row from next
        MPI_Irecv(grid[(rank + 1) * workload_rows + workload_rest - 1], GRID_SIZE, MPI_CXX_BOOL, next_proc, 0, MPI_COMM_WORLD, &req[0]);

        // send top row to previous process
        MPI_Isend(grid[rank * workload_rows], GRID_SIZE, MPI_CXX_BOOL, prev_proc, 0, MPI_COMM_WORLD, &req[1]);
        // printf("PROC %d, sending to %d\n", rank, prev_proc);

        // printf("PROC %d, recving to %d\n", rank, prev_proc);
        // recv on top row from previous process
        MPI_Irecv(grid[rank * workload_rows - 1], GRID_SIZE, MPI_CXX_BOOL, prev_proc, 0,
                  MPI_COMM_WORLD, &req[2]);

        // send bottom row to next process
        MPI_Isend(grid[(rank + 1) * workload_rows + workload_rest - 2], GRID_SIZE, MPI_CXX_BOOL, next_proc, 0, MPI_COMM_WORLD, &req[3]);

        MPI_Waitall(4, req, MPI_STATUS_IGNORE);

        // printf("PROC %d, sending to %d\n", rank, next_proc);
    }

    // Fix corners // MESSAGE PASSING
    // grid[0][0] = grid[GRID_SIZE - 2][GRID_SIZE - 2];
    // grid[GRID_SIZE - 1][GRID_SIZE - 1] = grid[1][1];
    // grid[0][GRID_SIZE - 1] = grid[GRID_SIZE - 2][1];
    // grid[GRID_SIZE - 1][0] = grid[1][GRID_SIZE - 2];

    if (rank == 0)
    {

        MPI_Request req[4];

        MPI_Irecv(&grid[0][0], GRID_SIZE, MPI_CXX_BOOL, size - 1, 0,
                  MPI_COMM_WORLD, &req[0]);

        MPI_Isend(&grid[1][1], 1, MPI_CXX_BOOL, size - 1, 0, MPI_COMM_WORLD, &req[1]);
        MPI_Irecv(&grid[0][GRID_SIZE - 1], GRID_SIZE, MPI_CXX_BOOL, size - 1, 0,
                  MPI_COMM_WORLD, &req[2]);
        MPI_Isend(&grid[1][GRID_SIZE - 2], 1, MPI_CXX_BOOL, size - 1, 0, MPI_COMM_WORLD, &req[3]);

        MPI_Waitall(4, req, MPI_STATUS_IGNORE);
    }
    if (rank == size - 1)
    {
        MPI_Request req[4];
        MPI_Isend(&grid[GRID_SIZE - 2][GRID_SIZE - 2], 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, &req[0]);

        MPI_Irecv(&grid[GRID_SIZE - 1][GRID_SIZE - 1], GRID_SIZE, MPI_CXX_BOOL, 0, 0,
                  MPI_COMM_WORLD, &req[1]);
        MPI_Isend(&grid[GRID_SIZE - 2][1], 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, &req[2]);
        MPI_Irecv(&grid[GRID_SIZE - 1][0], GRID_SIZE, MPI_CXX_BOOL, 0, 0,
                  MPI_COMM_WORLD, &req[3]);

        MPI_Waitall(4, req, MPI_STATUS_IGNORE);
    }

    // // columns padding must be copied anyhow (MESSAGE PASSING)
    // for (int j = 1; j < GRID_SIZE - 1; j++)
    // {
    //     // join columns together
    //     grid[0][j] = grid[GRID_SIZE - 2][j];
    //     grid[GRID_SIZE - 1][j] = grid[1][j];
    // }

    // printf("PROC %d Finished copyedeges\n", rank);
}

int readProblemFromInput(ProblemData &data)
{
    auto &grid = *data.readGrid;

    unsigned int seed = 0;
    std::cout << "READY" << std::endl;
    std::cin >> seed;

    std::cout << "Using seed " << seed << std::endl;
    if (seed == 0)
    {
        std::cout << "Warning: default value 0 used as seed." << std::endl;
    }

    // "random" numbers
    my_seedGenerator(seed);

    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i += 1)
    {
        *(grid[0] + i) = generateBit();
    }

    return seed;
}

void generateProblem(ProblemData &data, int seed)
{
    auto &grid = *data.readGrid;
    // printf("Process %d seed\n", seed);

    my_seedGenerator(seed);

    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i += 1)
    {
        *(grid[0] + i) = generateBit();
    }
}

int my_countAlive(ProblemData &data, int rank, int size)
{

    int workload_rows = GRID_SIZE / size;
    int workload_rest = GRID_SIZE % size;

    int flag = ((rank + 1) * workload_rows);

    if (rank == size - 1)
    {
        flag += workload_rest - 1;
    }

    auto &grid = *data.readGrid;
    int counter = 0;
    for (int x = rank == 0 ? 1 : rank * workload_rows; x < flag; x++)
    {
        for (int y = 1; y < GRID_SIZE - 1; y++)
        {
            if (grid[x][y])
            {
                counter++;
            }
        }
    }
    return counter;
}

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto *problemData = new ProblemData;

    int seed;
    if (rank == 0)
    {
        seed = readProblemFromInput(*problemData);
    }

    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // printf("Process %d seed\n", seed);
    if (rank != 0)
    {
        generateProblem(*problemData, seed);
    }

    // TODO@Students: This is the main simulation. Parallelize it using MPI.
    for (int iteration = 0; iteration < NUM_SIMULATION_STEPS; ++iteration)
    {

        {
            copy_edges(*problemData->readGrid, rank, size);

            if (iteration % SOLUTION_REPORT_INTERVAL == 0)
            {
                int res = 0;
                int local = my_countAlive(*problemData, rank, size);

                // printf("ITER %d, Local di %d is %d\n", iteration, rank, local);
                MPI_Reduce(&local, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if (rank == 0)
                {
                    std::cout << "Iteration " << iteration << ": " << res << " cells alive." << std::endl;
                }
            }

            evolve(*problemData, rank, size);

            problemData->swapGrids();
        }
    }

    int final_res = 0;
    int local = my_countAlive(*problemData, rank, size);
    // printf("ITER %d, Local di %d is %d\n", NUM_SIMULATION_STEPS, rank, local);
    MPI_Reduce(&local, &final_res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "Iteration " << NUM_SIMULATION_STEPS << ": " << final_res << " cells alive." << std::endl;
        std::cout << "DONE" << std::endl;
    }

    delete problemData;
    MPI_Finalize();
    return 0;
}
