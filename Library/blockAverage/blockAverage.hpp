#ifndef __blockAverage__
#define __blockAverage__

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "../PRNG/random.hpp"
#include "../library.hpp"
using namespace std;    

class BlockAverage{
    public:
        //constructor
        BlockAverage(int M_, int N_) : rnd("../../Library/PRNG/") {
            (M_ % N_ == 0) ? void() : (cerr << "ERROR: M must be a multiple of N" << endl, exit(1));
            M = M_;
            N = N_;
            L = M/N;
        }
        //destructor
        ~BlockAverage(){;}
        //methods
        void Run(std::ofstream& WriteData){
            double ave, ave2;
            double sum = 0, sum2 = 0;
            double sum_prog, sum2_prog;
            double err;

            for (int i = 0 ; i < N ; i++) {
                ave = block(); //computes a one-block average
                ave2 = ave*ave;

                sum += ave;
                sum2 += ave2;

                sum_prog = (double)sum/(i+1);
                sum2_prog = (double)sum2/(i+1);
                err = computeError(sum_prog, sum2_prog, i);
                WriteData << (i+1)*L << " " << sum_prog << " " << err << endl;
        }

        mean = sum_prog;
        meanError = err;
        }

        virtual double block() = 0;
        double GetMean() {return mean;}
        double GetError() {return meanError;}

        protected:
            int M, N, L;
            Random rnd; 
            double mean, meanError;

};

class Uniform_Mean : public BlockAverage{ 

    public:
        Uniform_Mean(int M, int N) : BlockAverage(M, N) {;}
        virtual ~Uniform_Mean() {;}
        double block() override{

            double val = 0;
            for (int i = 0 ; i < L ; i++) {val += rnd.Rannyu();}
            return val/L; 
        }
};

class Uniform_Sigma : public BlockAverage{ 

    public:
        Uniform_Sigma(int M, int N) : BlockAverage(M, N) {;}
        virtual ~Uniform_Sigma() {;}
        double block() override{

            double val = 0;
            for (int i=0 ; i<L ; i++) {val += pow(rnd.Rannyu()-0.5, 2);}
            return val/L; 
        }
};

class Buffon : public BlockAverage{ 

    public:
        Buffon(int M, int N, double L_, double d_) : BlockAverage(M, N), Length(L_), d(d_) {;}
        virtual ~Buffon() {;}
        double block() override{
            probability = 0;
            for (int j = 0; j < L; j++){
                double x = rnd.Rannyu(0, d/2.);
                double theta = rnd.Angle();
                if (x <= (Length/2)*fabs(sin(theta)))
                    probability += 1./L;
            }   
            return (2*Length)/(d*probability);
        }

    private:
        double Length, d, probability;
};

class IntegralUniform : public BlockAverage{ 

    public:
        IntegralUniform(int M, int N) : BlockAverage(M, N) {;}
        virtual ~IntegralUniform() {;}
        double block() override{
            double x = 0;
            for (int i = 0; i < L; i++){
            x += EvalCos(rnd.Rannyu())/L;
            }
            return x;
        }
};

class IntegralSampling : public BlockAverage{ 

    public:
        IntegralSampling(int M, int N) : BlockAverage(M, N) {;}
        virtual ~IntegralSampling() {;}
        double block() override{
            double x = 0;
            for (int i = 0; i < L; i++){
            x += EvalCos_sample(rnd.Sampling2_1())/L;
            }
            return x;
        }
};

/*class DiscreteWalker : public BlockAverage{ 

    public:
        DiscreteWalker(int M, int N) : BlockAverage(M, N) {;}
        virtual ~DiscreteWalker() {;}
        double block() override{
            int y;
            vector <double> x;
            for (int i = 0; i < L; i++){
            y = rnd.DiscreteMovement();
            x[y/2] += ((y % 2 == 0) ? 1. : -1.);
            }

        }
            
        
};*/


#endif
