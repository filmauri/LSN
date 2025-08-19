#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "../../../Library/PRNG/random.hpp"
#include "../../../Library/blockAverage/blockAverage.hpp"
#include "../../../Library/library.hpp"
#include "geneticLibrary.hpp"
#include "cities.hpp"

#include "mpi.h"

using namespace std;
int main(int argc, char* argv[]){

    // Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TSP tsp;
    cout << "Initializing TSP..." << endl;
    tsp.initialize(0);
    Random _rnd("../../../Library/PRNG/");

    arma::Mat <double> cities_positions;
    int nCities = tsp.getNCities();
    cities_positions.resize(nCities, 2);

    arma::Mat <int> migrators;
    migrators.resize(nCities, size);

    arma::Col <int> migrator;
    migrator.resize(nCities);

    if (rank == 0) 
        cities_positions = tsp.getCitiesPositions();

    MPI_Bcast(cities_positions.memptr(), nCities * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    cout << "Cities positions broadcasted." << endl;

    if (rank != 0) {
        tsp.setCitiesPositions(cities_positions);
        tsp.generate_tsp_population();
        tsp.sortPopulationbyLoss();
    }
    for (int i = 0; i < tsp.getNGenerations(); i++) {
        cout << "Generation " << i + 1 << endl;
        tsp.evolution(i);
        if (i % tsp.getMigrationStep() == 0 && i != 0) {
            cout << "Migrating individuals..." << endl;
            migrator = tsp.getbestTravel(i);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(migrator.memptr(), nCities, MPI_INT, migrators.memptr(), nCities, MPI_INT, 0, MPI_COMM_WORLD);
           
            if(rank==0)
                migrators = arma::shuffle(migrators);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(migrators.memptr(), nCities * 2, MPI_INT, 0, MPI_COMM_WORLD);

            tsp.set_best_travel(migrator);
            if (rank < 4)
                tsp.bestTravel(i);
           
            tsp.loss();
        }
        tsp.checkPopulation();
        tsp.checkStartingPos();
        tsp.saveBestLoss(i);
        if (i % 10 == 0 && rank == 0) {
            tsp.partialFinalization(i);
        }

        
    }
    tsp.finalize();
    MPI_Finalize();
    return 0;
}
