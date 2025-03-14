#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "../../Library/blockAverage/blockAverage.hpp"
#include "../../Library/PRNG/random.hpp"
#include "../../Library/library.hpp"

using namespace std;
int main(){
    Random rnd("../../Library/PRNG/");
    int M = 100000;              //Total number of throws
    int N = 100;                 // Number of blocks
    //int L = int(M/N);               //Number of throws in each block, please use for M a multiple of N
    Uniform_Mean test(M, N);
    ofstream WriteData;
    WriteData.open("dataClass.dat");
        test.Run(WriteData);
    WriteData.close();
    Uniform_Sigma test2(M, N);
    ofstream WriteSigma;
    WriteSigma.open("dataSigmaClass.dat");
        test2.Run(WriteSigma);
    WriteSigma.close();
    
    //chi quadro test
    int M2 = 100000; //Number of throws
    vector <double> chi2_values(N, 0); // Vector to store chi2 values

    for (int k = 0; k < N; k++) {
        vector<int> counts(N, 0); // Vector to count occurrences in each sub-interval

        // Generate n random numbers and count occurrences in each sub-interval
        for (int i = 0; i < M2; i++) {
            double r = rnd.Rannyu();
            int index = int(r * N);
            if (index < 0 || index >= N) {
                cerr << "ERROR: Index out of bounds: " << index << endl;
                return 1;
            }
            counts[index]++;
        }
        // Calculate chi2 for this set of random numbers
        double chi2 = 0;
        double expected = double(M2) / N;
        for (int i = 0; i < N; i++) {
            chi2 += pow(counts[i] - expected, 2) / expected;
        }
        chi2_values[k] = chi2;
    }
    ofstream WriteChi2;
    WriteChi2.open("chi2Class.dat");
    if (WriteChi2.is_open()) {
        for (int j = 0; j < N; j++) {
            WriteChi2 << j + 1 << " " << chi2_values[j] << endl;
        }
    }
   return 0;
}