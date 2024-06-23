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

    int flag = ((rank + 1) * workload_rows) - 1;

    if (rank == size - 1)
    {
        flag += workload_rest;
    }

    // For each cell
    for (int i = 1 + (rank * workload_rows); i < flag; i++)
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

    if (rank != size - 1)
    {
        // Copy data to the boundaries
        for (int i = 1 + (rank * workload_rows); i < ((rank + 1) * workload_rows) - 1; i++)
        {
            // join rows together
            grid[i][0] = grid[i][GRID_SIZE - 2];
            grid[i][GRID_SIZE - 1] = grid[i][1];
        }

        // send top row to previous process
        MPI_Send(grid + (rank * workload_rows), GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);

        // recv on top row from previous process
        MPI_Recv(grid + 1 + (rank * workload_rows), GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 0,
                 MPI_COMM_WORLD, nullptr);

        // send bottom row to next process
        MPI_Send(grid + ((rank + 1) * workload_rows), GRID_SIZE, MPI_CXX_BOOL, (rank + 1) % size, 1, MPI_COMM_WORLD);

        // recv on bottom row from next
        MPI_Recv(grid - 1 + ((rank + 1) * workload_rows), GRID_SIZE, MPI_CXX_BOOL, (rank + 1) % size, 1,
                 MPI_COMM_WORLD, nullptr);
    }
    else
    {
        printf("Last process of the line, should be 1\n");

        // Copy data to the boundaries
        for (int i = 1 + (rank * workload_rows); i < ((rank + 1) * workload_rows + workload_rest) - 1; i++)
        {
            // join rows together
            grid[i][0] = grid[i][GRID_SIZE - 2];
            grid[i][GRID_SIZE - 1] = grid[i][1];
        }
        // recv on top row from previous process
        MPI_Recv(grid + 1 + (rank * workload_rows), GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 0,
                 MPI_COMM_WORLD, nullptr);

        // send top row to previous process
        MPI_Send(grid + (rank * workload_rows), GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);

        // recv bottom padding row from next
        MPI_Recv(grid - 1 + ((rank + 1) * workload_rows) + workload_rest, GRID_SIZE, MPI_CXX_BOOL, (rank + 1) % size, 1, MPI_COMM_WORLD, nullptr);

        // send bottom row to next process
        MPI_Send(grid + ((rank + 1) * workload_rows) + workload_rest, GRID_SIZE, MPI_CXX_BOOL, (rank + 1) % size, 1, MPI_COMM_WORLD);
    }

    // Fix corners // MESSAGE PASSING
    // grid[0][0] = grid[GRID_SIZE - 2][GRID_SIZE - 2];
    // grid[GRID_SIZE - 1][GRID_SIZE - 1] = grid[1][1];
    // grid[0][GRID_SIZE - 1] = grid[GRID_SIZE - 2][1];
    // grid[GRID_SIZE - 1][0] = grid[1][GRID_SIZE - 2];

    if (rank == 0)
    {

        MPI_Recv(&grid[0][0], GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 0,
                 MPI_COMM_WORLD, nullptr);
        MPI_Recv(&grid[0][GRID_SIZE - 1], GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 1,
                 MPI_COMM_WORLD, nullptr);

        MPI_Send(&grid[1][1], 1, MPI_CXX_BOOL, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);
        MPI_Send(&grid[1][GRID_SIZE - 2], 1, MPI_CXX_BOOL, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);
    }
    else if (rank == size - 1)
    {
        MPI_Recv(&grid[GRID_SIZE - 1][GRID_SIZE - 1], GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 0,
                 MPI_COMM_WORLD, nullptr);
        MPI_Recv(&grid[GRID_SIZE - 1][0], GRID_SIZE, MPI_CXX_BOOL, (rank - 1 + size) % size, 1,
                 MPI_COMM_WORLD, nullptr);
        MPI_Send(&grid[GRID_SIZE - 2][GRID_SIZE - 2], 1, MPI_CXX_BOOL, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);
        MPI_Send(&grid[GRID_SIZE - 2][1], 1, MPI_CXX_BOOL, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);
    }

    // // columns padding must be copied anyhow (MESSAGE PASSING)
    // for (int j = 1; j < GRID_SIZE - 1; j++)
    // {
    //     // join columns together
    //     grid[0][j] = grid[GRID_SIZE - 2][j];
    //     grid[GRID_SIZE - 1][j] = grid[1][j];
    // }
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

    int flag = ((rank + 1) * workload_rows) - 1;

    if (rank == size - 1)
    {
        flag += workload_rest;
    }

    auto &grid = *data.readGrid;
    int counter = 0;
    for (int x = 1 + (rank * workload_rows); x < flag; x++)
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
    // printf("Process %d speaking\n", rank);
    if (rank == 0)
    {
        seed = readProblemFromInput(*problemData);
    }

    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Process %d seed\n", seed);
    if (rank == 0)
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
                MPI_Reduce(&local, &res, GRID_SIZE, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if (rank == 0)
                {
                    std::cout << "Iteration " << iteration << ": " << res << " cells alive." << std::endl;
                }
            }

            evolve(*problemData, rank, size);

            problemData->swapGrids();
        }
    }

    Utility::outputSolution(*problemData);

    delete problemData;
    MPI_Finalize();
    return 0;
}
