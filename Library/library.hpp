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
#endif