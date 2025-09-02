#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "../../Library/PRNG/random.hpp"
#include "../../Library/blockAverage/blockAverage.hpp"
#include "../../Library/library.hpp"

using namespace std;
int main()
{   
    double sigma = 1., mu = 2.;
    double x = 0.;
    int accepted = 1, attempted = 1;
    double integral = 0.;
    double step = 0.1;

    //cout << eval_H_psi(0.1, 0.0, 1) << endl;
    
    ofstream WriteData;
    WriteData.open("output.dat");
    ofstream WriteCoordinates;
    WriteCoordinates.open("coord.dat");
    int M = 100000, N = 100;
    int L = M / N;
    int nstep_equilibration = 100;
    Random rnd("../../Library/PRNG/");

    while ((double)accepted/attempted < 0.45 || (double)accepted/attempted > 0.55){
        step += 0.1;
        accepted = 0;
        attempted = 0;
        for (int i = 0; i < nstep_equilibration; i++){
            metro(x, rnd, sigma, mu, step, attempted, accepted);
        }
    }

    cout << "After the equilibration time, the step is: " << step << " with an acceptance of " << (double)accepted/attempted << endl;
    double ave = 0, ave2 = 0;
    double sum = 0, sum2 = 0;
    double sum_prog, sum2_prog;
    double err;

    for (int i = 0 ; i < N ; i++) {
        integral = 0;
        accepted = 0;
        attempted = 0;
        for (int j = 0; j < L; j++){
            metro(x, rnd, sigma, mu, step, attempted, accepted);
            integral += eval_H_psi(x, mu, sigma);
            WriteCoordinates << setw(8) << j + L * i
							 << setw(16) << x << endl;
            }
        ave = integral/(double)L; 
        ave2 = ave*ave;

        sum += ave;
        sum2 += ave2;

        sum_prog = sum/(i+1);
        sum2_prog = sum2/(i+1);
        
        err = computeError(sum_prog, sum2_prog, i);
        WriteData << (i + 1) * L << " " << sum_prog << " " << err << endl;
        }

    return 0;
}