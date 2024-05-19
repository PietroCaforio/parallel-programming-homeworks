//
// Created by Dennis-Florian Herr on 13/06/2022.
//

#include <string>
#include <deque>
#include <future>
#include <functional>

#include "Utility.h"
#include <iostream>

#define MEASURE_TIME true
#define NUM_WORKERS 32
#define NUM_GENERATORS 1

struct Problem {
    Sha1Hash sha1_hash;
    int problemNum;
};

std::atomic<int> finished_generator = 0;
std::atomic<int> work_done = 0;

int leadingZerosProblem = 8;
int leadingZerosSolution = 11;
int numProblems = 10000;
Sha1Hash *solutionHashes;
std::string * problemRandG;


/*
 * TODO@Students: Implement a thread safe queue.
 * Tip: use a condition variable to make threads wait when the queue is empty and there is nothing to pop().
 * https://en.cppreference.com/w/cpp/thread/condition_variable
 */
class ProblemQueue {
    public:
        void push(Problem problem){
            {
                std::lock_guard<std::mutex> lock(mutex);
                problem_queue.push_back(problem);
            }
            cv.notify_one();
        }

        Problem pop(){
            std :: unique_lock<std::mutex> lock ( mutex );
            while (problem_queue.empty()){
                cv.wait( lock );
            }
            Problem p = problem_queue.front();
            problem_queue.pop_front();
            return p;
        }

        bool empty(){
            return problem_queue.empty();
        }

    private:
        std::deque<Problem> problem_queue;
        std::mutex mutex;
        std::condition_variable cv;

};

ProblemQueue problemQueue;


// generate numProblems sha1 hashes with leadingZerosProblem leading zero bits
// This method is intentionally compute intense so you can already start working on solving
// problems while more problems are generated
void generateProblem(std::string problemRand[], int numProblems, int leadingZerosProblem, int idx){
    // std::cout << numProblems << "\n";
    Sha1Hash mock;

    int step = numProblems / NUM_GENERATORS;
    if (! (numProblems % NUM_GENERATORS == 0)) {std::cout << "ass error";exit(1);}

    for(int i = idx*step; i < (idx+1)*step; i++){
        
        Sha1Hash hash = Utility::sha1(problemRand[i]);
        do{
            // we keep hashing ourself until we find the desired amount of leading zeros
            hash = Utility::sha1(hash);
        }while(Utility::count_leading_zero_bits(hash) < leadingZerosProblem);
        problemQueue.push(Problem{hash, i});
    }

    finished_generator++;
    // std::cout << "GENERATOR FINISHED\n"; 
    if (finished_generator >= NUM_GENERATORS) {
        for (int i = 0; i < NUM_WORKERS; i++) {
            problemQueue.push(Problem{mock, -10});
        }
    }
}

// This method repeatedly hashes itself until the required amount of leading zero bits is found
Sha1Hash findSolutionHash(Sha1Hash hash, int leadingZerosSolution){
    do{
        // we keep hashing ourself until we find the desired amount of leading zeros
        hash = Utility::sha1(hash);
    }while(Utility::count_leading_zero_bits(hash) < leadingZerosSolution);

    return hash;
}

void worker() {
    while(true) {
        Problem p = problemQueue.pop();
        // std::cout << p.problemNum << "\n";
        if (p.problemNum < 0) {
            break;
        }
        solutionHashes[p.problemNum] = findSolutionHash(p.sha1_hash, leadingZerosSolution);
        work_done++;
    }
}
    

int main(int argc, char *argv[]) {
    

    //Not interesting for parallelization
    Utility::parse_input(numProblems, leadingZerosProblem, leadingZerosSolution, argc, argv);
    Sha1Hash solutionHashes2[numProblems];
    solutionHashes = new Sha1Hash[numProblems];
    problemRandG = new std::string[numProblems];
    
    unsigned int seed = Utility::readInput();

    #if MEASURE_TIME
    struct timespec generation_start, generation_end;
    clock_gettime(CLOCK_MONOTONIC, &generation_start);
    #endif

    /*
    * TODO@Students: Generate the problem in another thread and start already working on solving the problems while the generation continues
    */
    srand(seed+1);
    for (int i = 0; i < numProblems; i++) {
        problemRandG[i] = std::to_string(rand()) + std::to_string(rand());
    }

    std::thread generators [NUM_GENERATORS];
    for (int i = 0; i < NUM_GENERATORS; i++) {
        generators[i] = std::thread(generateProblem, problemRandG, numProblems, leadingZerosProblem, i);
    }

    #if MEASURE_TIME
    clock_gettime(CLOCK_MONOTONIC, &generation_end);
    double generation_time = (((double) generation_end.tv_sec + 1.0e-9 * generation_end.tv_nsec) - ((double) generation_start.tv_sec + 1.0e-9 * generation_start.tv_nsec));
    fprintf(stderr, "Generate Problem time:  %.7gs\n", generation_time);

    struct timespec solve_start, solve_end;
    clock_gettime(CLOCK_MONOTONIC, &solve_start);
    #endif

    /*
    * TODO@Students: Create worker threads that parallelize this functionality. Add the synchronization directly to the queue
    */
    //worker();
    std::thread workers[NUM_WORKERS];
    for (int i=0; i < NUM_WORKERS; i++) {
        workers[i] = std::thread(worker);
    }

    for (int i = 0; i < NUM_GENERATORS; i++) {
        generators[i].join();
    }
    
    for (int i = 0; i < NUM_WORKERS; i++) {
        workers[i].join();
    }
    
    #if MEASURE_TIME
    clock_gettime(CLOCK_MONOTONIC, &solve_end);
    double solve_time = (((double) solve_end.tv_sec + 1.0e-9 * solve_end.tv_nsec) - ((double) solve_start.tv_sec + 1.0e-9 * solve_start.tv_nsec));
    fprintf(stderr, "Solve Problem time:     %.7gs\n", solve_time);
    #endif

    /*
    * TODO@Students: Make sure all work has finished before calculating the solution
    * Tip: Push a special problem for each thread onto the queue that tells a thread to break and stop working
    */

    Sha1Hash solution;
    // guarantee initial solution hash data is zero
    memset(solution.data, 0,  SHA1_SIZE);
    // this doesn't need parallelization. it's neglectibly fast
    for(int i = 0; i < numProblems; i++){
        solution = Utility::sha1(solution, solutionHashes[i]);
    }

    Utility::printHash(solution);
    printf("DONE\n");

    return 0;
}
