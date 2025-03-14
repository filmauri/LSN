#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "../../Library/PRNG/random.hpp"
using namespace std;

double error(vector <double> av, vector <double> av2, int n)
{
    if (n == 0)
        return 0;
    else
        return sqrt((av2[n] - pow(av[n],2)) / n);
}

int main()
{
    Random rnd("../../Library/PRNG/");

    double M = 10000;              //Total number of throws
    double N = 100;                 //# Number of blocks
    int L = int(M/N);               //Number of throws in each block, please use for M a multiple of N
    //np.random.seed(1)             //Fixing random seed for reproducibility
    vector <double> x;              //[0,1,2,...,N-1]
    for(int i = 0; i < N; i++)
        x.push_back(i);

    vector <double> ave(N,0), aveSigma(N,0);;       //mean
    vector <double> av2(N,0), av2Sigma(N,0);       //mean of the squares
    vector <double> sum_prog(N,0), sum_progSigma(N,0);  //cumulative mean
    vector <double> su2_prog(N,0), su2_progSigma(N,0); //cumulative mean of the squares
    vector <double> err_prog(N,0), err_progSigma(N,0);  //statistical uncertainty

    double somma_prog = 0, somma_progSigma = 0; //Cumulative sum
    double somma2_prog = 0, somma2_progSigma = 0; //Cumulative sum of the squares
    double sum = 0, sumSigma = 0; //Sum of the throws

    for(int i = 0; i < N; i++)
    {   
        sum = 0;
        sumSigma = 0;
        for(int j = 0; j < L; j++){
            double y = rnd.Rannyu();
            sum += y;
            sumSigma += pow(y - 0.5, 2);
        }
        
        //mean value of the integral    
        ave[i] = sum/L; //Data blocking, mean of each block
        av2[i] = pow(ave[i],2);

        somma_prog += ave[i]; 
        somma2_prog += av2[i]; 
        sum_prog[i] = somma_prog/(i+1); //Cumulative mean
        su2_prog[i] = somma2_prog/(i+1); //Cumulative mean of the squares
        err_prog[i] = error(sum_prog, su2_prog, i); //Statistical uncertainty

        //mean value of deviation from the mean
        aveSigma[i] = sumSigma/L; 
        av2Sigma[i] = pow(aveSigma[i],2);

        somma_progSigma += aveSigma[i]; 
        somma2_progSigma += av2Sigma[i]; 
        sum_progSigma[i] = somma_progSigma/(i+1); 
        su2_progSigma[i] = somma2_progSigma/(i+1); 
        err_progSigma[i] = error(sum_progSigma, su2_progSigma, i); 

    }


    ofstream WriteData;
    WriteData.open("data.dat");
    if (WriteData.is_open()){
        for(int i = 0; i < N; i++)
            WriteData << x[i] * L << " " << sum_prog[i] - 0.5 << " " << err_prog[i] << endl;
    } else cerr << "PROBLEM: Unable to create data.dat" << endl;

    ofstream WriteSigma;
    WriteSigma.open("dataSigma.dat");
    if (WriteSigma.is_open()){
        for(int i = 0; i < N; i++)
            WriteSigma << x[i] * L << " " << sum_progSigma[i] - 1./12. << " " << err_progSigma[i] << endl;
    } else cerr << "PROBLEM: Unable to create dataSigma.dat" << endl;
    
    
    //chi quadro test
    int M2 = 100; //Number of sub-intervals
    int n = 100000; //Number of throws
    vector <double> chi2_values(M2, 0); // Vector to store chi2 values

    for (int k = 0; k < 100; k++) {
        vector<int> counts(M2, 0); // Vector to count occurrences in each sub-interval

        // Generate n random numbers and count occurrences in each sub-interval
        for (int i = 0; i < n; i++) {
            double r = rnd.Rannyu();
            int index = int(r * M2);
            if (index < 0 || index >= M2) {
                cerr << "ERROR: Index out of bounds: " << index << endl;
                return 1;
            }
            counts[index]++;
        }

        // Calculate chi2 for this set of random numbers
        double chi2 = 0;
        double expected = double(n) / M2;
        for (int i = 0; i < M2; i++) {
            chi2 += pow(counts[i] - expected, 2) / expected;
        }
        chi2_values[k] = chi2;
    }

    ofstream WriteChi2;
    WriteChi2.open("chi2.dat");
    if (WriteChi2.is_open()) {
        for (int j = 0; j < M2; j++) {
            WriteChi2 << j + 1 << " " << chi2_values[j] << endl;
        }
    } else {
        cerr << "PROBLEM: Unable to create chi2.dat" << endl;
    }

    
    return 0;
}