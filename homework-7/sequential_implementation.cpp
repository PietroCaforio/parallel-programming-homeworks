//
// Created by Vincent Bode on 08/07/2020.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "life.h"
#include "Utility.h"
#include "VideoOutput.h"

/*
  Apply the game of life rules on a Torus --> grid contains shadow rows and columns
  to simplify application of rules i.e. grid actually ranges from grid [ 1.. height - 2 ][ 1 .. width - 2]
*/
void evolve(ProblemData &problemData)
{
    auto &grid = *problemData.readGrid;
    auto &writeGrid = *problemData.writeGrid;
    // For each cell
    for (int i = 1; i < GRID_SIZE - 1; i++)
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
void copy_edges(bool (&grid)[GRID_SIZE][GRID_SIZE])
{
    // Copy data to the boundaries
    for (int i = 1; i < GRID_SIZE - 1; i++)
    {
        // join rows together
        grid[i][0] = grid[i][GRID_SIZE - 2];
        grid[i][GRID_SIZE - 1] = grid[i][1];
    }

    for (int j = 1; j < GRID_SIZE - 1; j++)
    {
        // join columns together
        grid[0][j] = grid[GRID_SIZE - 2][j];
        grid[GRID_SIZE - 1][j] = grid[1][j];
    }

    // Fix corners
    grid[0][0] = grid[GRID_SIZE - 2][GRID_SIZE - 2];
    grid[GRID_SIZE - 1][GRID_SIZE - 1] = grid[1][1];
    grid[0][GRID_SIZE - 1] = grid[GRID_SIZE - 2][1];
    grid[GRID_SIZE - 1][0] = grid[1][GRID_SIZE - 2];
}

int my_countAlive(ProblemData &data, int rank, int size)
{

    int workload_rows = GRID_SIZE / size;
    int workload_rest = GRID_SIZE % size;

    printf("Workload %d, rest %d\n", workload_rows, workload_rest);

    int flag = ((rank + 1) * workload_rows);

    if (rank == size - 1)
    {
        flag += workload_rest - 1;
    }

    printf("Tarn %d,  Flag %d \n", rank == 0 ? 1 : rank * workload_rows, flag);

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

    auto *problemData = new ProblemData;

    // As with Jack Sparrow's exercise, this needs FFMPEG (new and improved: this now works with more video players).
    // As an alternative, you can write individual png files to take a look at the data.

    Utility::readProblemFromInput(*problemData);

    // TODO@Students: This is the main simulation. Parallelize it using MPI.
    for (int iteration = 0; iteration < NUM_SIMULATION_STEPS; ++iteration)
    {
        if (iteration % SOLUTION_REPORT_INTERVAL == 0)
        {
            printf("MC live %d --- %d  SUM %d\n", my_countAlive(*problemData, 0, 2), my_countAlive(*problemData, 1, 2), my_countAlive(*problemData, 0, 2) + my_countAlive(*problemData, 1, 2));
            Utility::outputIntermediateSolution(iteration, *problemData);
        }

        copy_edges(*problemData->readGrid);

        if (iteration % SOLUTION_REPORT_INTERVAL == 0)
        {
            Utility::outputIntermediateSolution(iteration, *problemData);
        }

        evolve(*problemData);

        problemData->swapGrids();
    }

    Utility::outputSolution(*problemData);

    delete problemData;
    return 0;
}
