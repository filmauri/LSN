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

    double energy = 0.;
    double energy_Old = energy;

    double mu = 1., sigma = 1.;
    double mu_Old = mu, sigma_Old = sigma;

    double x = 0.;
    double integral = 0.;
    double step = 0.005;

    // File per tutti i grafici
    ofstream WriteData;
    WriteData.open("SA_history.dat");  // File principale per SA history
    ofstream WriteParameters;
    WriteParameters.open("parameters_trajectory.dat");  // Traiettoria parametri
    
    int M = 100000, N = 100;
    int L = M / N;
    int accepted = 1, attempted = 1;
    int SA_step = 0;  // Contatore passi SA
    
    // Headers dei file
    WriteData << "# SA_step T energy error mu sigma acceptance" << endl;
    WriteParameters << "# SA_step mu sigma T energy" << endl;
    
    x = 0.0;
    accepted = 1; 
    attempted = 1;

    cout << "Inizio Simulated Annealing..." << endl;
    beta = 0.1;
    while(SA_step < 1000){
        beta += 0.1;
        T = 1./beta;
        SA_step++;
        
        // Proponi nuovi parametri
        double mu_new = abs(mu + rnd.Gauss(0,step*T));
        double sigma_new = abs(sigma + rnd.Gauss(0,step*T));
        
        while(sigma_new < 0.2) {
            sigma_new = abs(sigma + rnd.Rannyu() * T);
        }

        // Calcola energia con i nuovi parametri
        integral = 0;
        double ave = 0, ave2 = 0;
        double sum = 0, sum2 = 0;
        double sum_prog, sum2_prog;
        double err;
        double para_acceptance = 0;
        double metrostep = 0.05;  // Mantieni step fisso
        
        // Reset posizione per ogni calcolo energia
        x = 0.0;
        accepted = 0;
        attempted = 0;
        
        for (int i = 0 ; i < N ; i++) {
            integral = 0;
            for (int j = 0; j < L; j++){
                metro(x, rnd, sigma_new, mu_new, metrostep, attempted, accepted);
                integral += eval_H_psi(x, mu_new, sigma_new);
            }
            ave = integral/(double)L; 
            ave2 = ave * ave;

            sum += ave;
            sum2 += ave2;

            sum_prog = sum/(i+1);
            sum2_prog = sum2/(i+1);
        
            err = computeError(sum_prog, sum2_prog, i);
        }
        
        energy = sum_prog;
        
        para_acceptance = min(1., exp(-beta * (energy - energy_Old)));
        
        bool accepted_step = false;
        if (rnd.Rannyu() < para_acceptance){
            // Accetta: aggiorna tutto
            energy_Old = energy;
            mu = mu_new;
            sigma = sigma_new;
            accepted_step = true;
        }
        
        // Scrivi nei file (sempre i parametri ATTUALI)
        WriteData   << setw(8) << SA_step
                    << setw(15) << T 
                    << setw(15) << energy_Old  // Energia accettata
                    << setw(15) << err
                    << setw(15) << mu 
                    << setw(15) << sigma 
                    << setw(15) << (double)accepted/attempted << endl;
                    
        WriteParameters << setw(8) << SA_step
                       << setw(15) << mu
                       << setw(15) << sigma  
                       << setw(15) << T
                       << setw(15) << energy_Old << endl;
    }
    
    WriteData.close();
    WriteParameters.close();
    
    cout << "\nSimulated Annealing completato!" << endl;
    cout << "Parametri finali:" << endl;
    cout << "T finale: " << T << endl;
    cout << "Mu finale: " << mu << endl;
    cout << "Sigma finale: " << sigma << endl;
    cout << "Energia finale: " << energy_Old << endl;

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
    
    int N_final = 200;  
    int L_final = M / N_final;
    
    x = 0.0;
    accepted = 0;
    attempted = 0;
    double sum_final = 0, sum2_final = 0, ave = 0;
    
    for (int i = 0; i < N_final; i++) {
        
        integral = 0;
        for (int j = 0; j < L_final; j++){
            metro(x, rnd, sigma, mu, step, attempted, accepted);
            
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