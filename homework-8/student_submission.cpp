#include <cstddef> // for size_t
#include "Utility.h"
#include <algorithm>
#include <mpi.h>

double local_min{0.0};
double local_max{0.0};

void compute_stencil(ProblemData &pd)
{
    // For evert cell with coordinates (y,x) compute the influx from neighbor cells
    // Apply reflecting boundary conditions
    double **domain = Utility::get_domain(pd);

    for (size_t y = 0; y < pd.dimension; y++)
    {
        for (size_t x = 0; x < pd.dimension; x++)
        {

            double cell_water = domain[y][x];
            double update = 0.0;

            // Add left neighbor
            if (x != 0)
            {
                double difference = domain[y][x - 1] - cell_water;
                update += difference / Utility::viscosity_factor;
            }

            // Add right neighbor
            if (x != pd.dimension - 1)
            {
                double difference = domain[y][x + 1] - cell_water;
                update += difference / Utility::viscosity_factor;
            }

            // Add lower neighbor
            if (y != 0)
            {
                double difference = domain[y - 1][x] - cell_water;
                update += difference / Utility::viscosity_factor;
            }

            // Add upper neighbor
            if (y != pd.dimension - 1)
            {
                double difference = domain[y + 1][x] - cell_water;
                update += difference / Utility::viscosity_factor;
            }

            double waterheight = domain[y][x] + update;

            if (pd.use_second_domain)
            {
                pd.domain1[y][x] = waterheight;
            }
            else
            {
                pd.domain2[y][x] = waterheight;
            }
        }
    }
}

bool termination_criteria_fulfilled(ProblemData &pd)
{
    double max = local_max;
    double min = local_min;

    // TODO @Students:
    // track min/max across ranks

    if ((max - min) < Utility::threshold)
    {
        return true;
    }

    return false;
}

void exchange_halo(ProblemData &pd)
{
    // TODO @Students:
    // Implement halo exchange
}

unsigned long long simulate(ProblemData &pd)
{
    volatile bool terminate_criteria_met = false;

    while (!terminate_criteria_met)
    {
        exchange_halo(pd);
        compute_stencil(pd);
        terminate_criteria_met = termination_criteria_fulfilled(pd);
        Utility::switch_arrays(pd);
        pd.patch_updates += 1;
    }

    return pd.patch_updates;
}

double **generate_initial_water_height(ProblemData &pd, int seed)
{
    Utility::generator.seed(seed);
    size_t half_dimension = pd.dimension / 2;
    size_t y_offsets[2] = {0, half_dimension};
    size_t x_offsets[2] = {0, half_dimension};

    constexpr size_t total_domain_count = 4;
    double **data = new double *[total_domain_count];

    for (size_t i = 0; i < total_domain_count; i++)
    {
        data[i] = new double[pd.dimension * pd.dimension];
    }

    size_t domain_counter = 0;
    for (size_t yoff : y_offsets)
    {
        for (size_t xoff : x_offsets)
        {
            for (size_t y = 0 + yoff; y < half_dimension + yoff; y++)
            {
                for (size_t x = 0 + xoff; x < half_dimension + xoff; x++)
                {
                    data[domain_counter][y * pd.dimension + x] = Utility::get_water_height(y, x);
                }
            }

            domain_counter += 1;
        }
    }

    return data;
}

int summer(ProblemData &pd, int size) {
    int sum = 0;
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++){
            
            sum += pd.domain1[i][j];
            
        }
    }
    return sum;
}


int main(int argc, char* argv[])
{
    ProblemData pd;

    // TODO @Students:
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int seed = 0;
    if (my_rank == 0) {
        std::cout << "READY" << std::endl;
        std::cin >> seed;    
    }
    int my_size = (Utility::domain_size / size );
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    
    // TODO @Students:
    // Think about minimizing the array size at each rank. 
    // this might require additional changes elsewhere
    
    Utility::init_problem_data(pd, my_size);
    double **initial_water_heights;
    if (my_rank == 0) {
        **initial_water_heights = generate_initial_water_height(pd, seed);
    }
    MPI_Scatter();
    int my_sum = summer(pd, my_size);
    int res = 0;
    
    
    // TODO @Students:
    // Initialize MPI find a way to send the initial data to the domains of other MPI Ranks
    double *this_rank_water_heights;
    Utility::apply_initial_water_height(pd, this_rank_water_heights);
    MPI_Reduce(&my_sum, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    printf("%d\n",res);
    MPI_Finalize();
    exit(0);
    delete[] initial_water_heights;

    simulate(pd);

    Utility::free_problem_data(pd);

    std::cout << pd.patch_updates << std::endl;

    // TODO @Students:
    // Finalize MPI
    

    return 0;
}
