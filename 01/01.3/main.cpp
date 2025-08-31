#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "../../Library/PRNG/random.hpp"
#include "../../Library/library.hpp"
#include "../../Library/blockAverage/blockAverage.hpp"
using namespace std;
int main(){
    Random rnd("../../Library/PRNG/");

    //Buffon's experiment
    double M = 1000000; //Number of throws
    double N = 100; //Number of blocks
    double L = 1.; //Needle length
    double d = 2.; //Distance between lines
    
    Buffon test(M, N, L, d);
    ofstream WriteBuffon;
    WriteBuffon.open("dataBuffonClass.dat");
    test.Run(WriteBuffon);
    WriteBuffon.close();
    return 0;
}