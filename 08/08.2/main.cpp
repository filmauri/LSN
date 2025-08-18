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
    Random rnd("../../Library/PRNG/");

    double T = 1; 
    double T_end = 0.001;
    double beta = 1./T;

    double energy = 1.;
    double energy_Old = energy;

    double mu = 0., sigma = 1.;
    double mu_Old = mu, sigma_Old = sigma;

    double x = 0.;
    double integral = 0.;
    double step = 0.1;

    ofstream WriteData;
    WriteData.open("output.dat");
    int M = 100000, N = 100;
    int L = M / N;
    int nstep_equilibration = 100;
    int accepted = 1, attempted = 1;

    while(T > T_end){
        T *= 0.995;
        beta = 1./T;
        mu = abs(mu + rnd.Rannyu(-1,1) * T);
        sigma = abs(sigma + rnd.Rannyu(-1,1) * T);

        integral = 0;
        double ave = 0, ave2 = 0;
        double sum = 0, sum2 = 0;
        double sum_prog, sum2_prog;
        double err;
        double para_acceptance = 0;
        step = 0.1;
        //equilibration and step setting 
        /*while ((double)accepted/attempted < 0.45 || (double)accepted/attempted > 0.55){
            step += 0.1;
            accepted = 0;
            attempted = 0;
            for (int i = 0; i < nstep_equilibration; i++){
                metro(x, rnd, sigma, mu, step, attempted, accepted);
            }
        }
        cout << "step " << step << endl; */
        
        for (int i = 0 ; i < N ; i++) {
            integral = 0;
            accepted = 0;
            attempted = 0;
            for (int j = 0; j < L; j++){
                metro(x, rnd, sigma, mu, step, attempted, accepted);
                integral += eval_H_psi(x, mu, sigma);
            }
            ave = integral/(double)L; 
            ave2 = ave * ave;

            sum += ave;
            sum2 += ave2;

            sum_prog = sum/(i+1);
            sum2_prog = sum2/(i+1);
        
            err = computeError(sum_prog, sum2_prog, i);
            energy = sum_prog;
        
            para_acceptance = min(1., exp(-beta * (energy - energy_Old)));
            
            if (rnd.Rannyu() < para_acceptance){
			energy_Old = energy;
			mu_Old = mu; 
			sigma_Old = sigma;
            }

            WriteData   << setw(15) << T 
                        << setw(15) << (double)accepted/attempted
                        << setw(15) << mu 
                        << setw(15) << sigma 
                        << setw(15) << energy 
                        << setw(15) << err << endl;
        }
    }
    cout << "T " << T << "\nMu "<< mu << "\nSigma " << sigma << endl;

    return 0;
}