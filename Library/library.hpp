#ifndef __LIBRARY__
#define __LIBRARY__
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
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
#endif