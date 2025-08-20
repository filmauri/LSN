#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "../../../../Library/PRNG/random.hpp"
#include "../../../../Library/library.hpp"
#include "geneticLibrary.hpp"
#include "cities.hpp"

using namespace std;
int main(){
        TSP tsp;
        cout << "Initializing TSP..." << endl;
        tsp.initialize();
        for (int i = 0; i < tsp.getNGenerations(); i++) {
            cout << "Generation " << i + 1 << endl;
            tsp.evolution(i);
            //tsp.printPopulation();
            //tsp.bestTravel(i);
            tsp.averageLoss(i);

            if (i % 10 == 0) 
                tsp.partialFinalization(i);
            
        }
        tsp.finalize();
    return 0;
}