#ifndef __LIBRARY__
#define __LIBRARY__
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

double computeError(double av, double av2, int n)
{
    if (n == 0)
        return 0;
    else
        return sqrt((av2 - pow(av,2)) / n);
}
double EvalCos(double x){return M_PI/2*cos(M_PI*x/2);}
double EvalCos_sample(double x){return M_PI/(4*(1-x))*cos(M_PI*x/2);}

double S_tot(double t_tot, double S0, double r, double sigma, double W){
    return S0 * exp((r - 0.5 * sigma * sigma) * t_tot + sigma * W);
}
double S_step(double t_step, double r, double sigma, double W){
    return exp((r - 0.5 * sigma * sigma) * t_step + sigma * W * sqrt(t_step));
}

double call_price(double r, double t_tot, double S, double K){
    return exp(-r * t_tot) * max(S - K, 0.0);
}

double put_price(double r, double t_tot, double S, double K){
    return exp(-r * t_tot) * max(K - S, 0.0);
}

//functions for lesson 8
double eval_Psi_2(double x, double sigma, double mu)
{
    double psi = exp(-pow((x - mu), 2)/(2 * pow(sigma, 2))) + exp(-pow((x + mu), 2)/(2 * pow(sigma, 2)));
    return psi * psi;
}

void metro(double &x, Random &rnd, double sigma, double mu, double step, int &attempted, int &accepted)
{
    double _x_new;
    _x_new = x + rnd.Rannyu(-1,1) * step;
    double alpha = min(1., eval_Psi_2(_x_new, sigma, mu)/eval_Psi_2(x, sigma, mu));
    if(rnd.Rannyu() < alpha){
        x = _x_new;
        accepted++;
    }
    attempted++;
}

double eval_H_psi(double x, double mu, double sigma) {
    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;

    double exp1 = exp(- (x - mu) * (x - mu) / (2.0 * sigma2));
    double exp2 = exp(- (x + mu) * (x + mu) / (2.0 * sigma2));
    double psi = exp1 + exp2;

    double d2_exp1 = (((x - mu)*(x - mu)) / sigma4) - (1./sigma2);
    double d2_exp2 = (((x + mu)*(x + mu)) / sigma4) - (1./sigma2);
    
    double der = exp1 * d2_exp1 + exp2 * d2_exp2;

    double kinetic = - 0.5 * der / psi;
    double potential = pow(x, 4) - 5./2. * x * x;

    return kinetic + potential;
}


#endif