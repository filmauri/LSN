#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "../../Library/PRNG/random.hpp"
#include "../../Library/blockAverage/blockAverage.hpp"
#include "../../Library/library.hpp"

using namespace std;
int main()
{   
    Random rnd("../../Library/PRNG/");

    double SA_step = 1000;
    double T = 1; 
    double beta = 1;
    double counter = 0.;

    double energy = 0.;
    double energy_Old = 10;
    double best_energy = energy;
    double energyError = 0.;
    double best_energyError = 0.;

    double mu = 1., sigma = 1.;
    double mu_Old = mu, sigma_Old = sigma;
    double bestMu, bestSigma;
    double min_step = 1;
    double metrostep = 3;
    double M = 100000;
    double N = 100;
    int L = M/N;
    double deltaMu = 0, deltaSigma = 0;
    int accepted=0, attempted=0;
    double x = 0.;


    double integral = 0;
    double ave = 0, ave2 = 0;
    double sum = 0, sum2 = 0;
    double sum_prog, sum2_prog;


    //file per primo grafico 
    ofstream fout("history_H.dat");
    fout << "#      step        H        error" << endl;
    //file per spazio parametri
    ofstream foutParameters("MuSigma.dat");
    foutParameters << "#      step        mu        sigma" << endl;

    //SA
    for (int i = 0; i < SA_step; i++)
    {   
        beta +=0.1;
        T = 1./beta;

        deltaMu = deltaMu > mu_Old ? mu_Old : min_step;
        deltaSigma = deltaSigma > sigma_Old ? sigma_Old : min_step;

        accepted = 0;
        attempted = 0;
        ave = 0, ave2 = 0;
        sum = 0, sum2 = 0;
        sum_prog = 0, sum2_prog = 0;
        //candido i nuovi parametri
        mu = mu_Old + rnd.Gauss(0, deltaMu * T);
        sigma = sigma_Old + rnd.Gauss(0, deltaSigma * T);

        //stimo il valore di H
        for (int index = 0 ; index < N ; index++) {
            integral = 0;
            for (int j = 0; j < L; j++){
                metro(x, rnd, sigma, mu, metrostep, attempted, accepted);
                integral += eval_H_psi(x, mu, sigma);
            }
            ave = integral/(double)L; 
            ave2 = ave * ave;

            sum += ave;
            sum2 += ave2;

            sum_prog = sum/(index+1);
            sum2_prog = sum2/(index+1);
        
            energyError = computeError(sum_prog, sum2_prog, index);
        }
        energy = sum_prog;

        //cout << "Accettanza " << (double)accepted/attempted << endl;     
        fout << i << setw(20) << energy << setw(20) << energyError << endl;
        foutParameters << i << setw(20) << mu << setw(20) << sigma << endl;
        
        //accettazione mossa parametri
        double alpha = min(1., exp(-beta * (energy - energy_Old)));
        if(rnd.Rannyu() < alpha){
            cout << "ciao" << endl;
            counter++;
            mu_Old = mu;
            sigma_Old = sigma;
            energy_Old = energy;
            if (energy < best_energy){   
                best_energy = energy;
                best_energyError = energyError;
                bestMu = mu;
                bestSigma = sigma;
                }
                
        }else{
            mu = mu_Old;
            sigma = sigma_Old;
            
        }   
    }
    cout << (double)counter/SA_step << endl;

    cout << "\nSimulated Annealing completato!" << endl;
    cout << "Parametri finali:" << endl;
    cout << "T finale: " << T << endl;
    cout << "Mu finale: " << bestMu << endl;
    cout << "Sigma finale: " << bestSigma << endl;
    cout << "Energia finale: " << best_energy << endl;
    mu = bestMu;
    sigma=bestSigma;

    // ============================================
    // SIMULAZIONE FINALE CON PARAMETRI OTTIMALI
    // ============================================
    
    cout << "\nInizio simulazione finale con parametri ottimali..." << endl;
    
    x = 0.0;
    accepted = 0;
    attempted = 0;
    
    ofstream WriteConvergence;
    WriteConvergence.open("convergence.dat");
    WriteConvergence << "# block_number total_steps progressive_average error" << endl;
    
    ofstream WriteSamples;
    WriteSamples.open("sampled_positions.dat");
    WriteSamples << "# x_positions" << endl;
    int M_final = 1000000;
    int N_final = 200;  
    int L_final = M_final / N_final;
    
    x = 0.0;
    accepted = 0;
    attempted = 0;
    ave =0;
    double sum_final = 0, sum2_final = 0;
    
    for (int i = 0; i < N_final; i++) {
        
        integral = 0;
        for (int j = 0; j < L_final; j++){
            metro(x, rnd, sigma, mu, metrostep, attempted, accepted);
            
            // Salva posizione per istogramma (ogni 10 punti per non avere file troppo grossi)
            if(j % 10 == 0) {
                WriteSamples << setw(15) << x << endl;
            }
            
            integral += eval_H_psi(x, mu, sigma);
        }
        
        ave = integral / (double)L_final;
        sum_final += ave;
        sum2_final += ave * ave;
        
        // Media progressiva
        double prog_ave = sum_final / (i + 1);
        double prog_av2 = sum2_final / (i + 1);
        double prog_err = computeError(prog_ave, prog_av2, i);

        WriteConvergence << setw(8) << (i+1) 
                        << setw(12) << (i+1) * L_final
                        << setw(15) << prog_ave
                        << setw(15) << prog_err << endl;
    }
    
    WriteConvergence.close();
    WriteSamples.close();
    
    
    return 0; 
}
