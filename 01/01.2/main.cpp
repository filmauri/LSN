#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "../../Library/PRNG/random.hpp"
using namespace std;
int main()
{   
    Random rnd("../../Library/PRNG/"); //inizializzazione di un elemento della classe Random
    // Creazione di un vettore di ofstream di 16 elementi
    vector<ofstream> output_files(4);
    for (int i = 0; i < 4; i++) {
        string filename = "output_" + to_string(i) + ".dat";
        output_files[i].open(filename);
        if (!output_files[i].is_open()) {
            cerr << "PROBLEM: Unable to open " << filename << endl;
        }
    }
    vector <int> N = {1, 2, 10, 100};
    int count = 0;

    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 10000; k++) {
            double meanU = 0;
            double meanExp = 0;
            double meanLor = 0;
            for (int j = 0; j < N[count]; j++) {
                meanU += rnd.Rannyu()/N[count];
                meanExp += rnd.Exponential(1)/N[count];
                meanLor += rnd.CauchyLorentz(0, 1)/N[count];
            }
            output_files[i] << meanU << " " << meanExp << " " << meanLor << endl;
        }
        count ++;
    }
    return 0;
}