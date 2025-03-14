#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "../../Library/PRNG/random.hpp"
#include "../../Library/blockAverage/blockAverage.hpp"
#include "../../Library/library.hpp"

using namespace std;
int main()
{   
    Random rnd("../../Library/PRNG/");
    int M = 100000; // Number of throws
    int N = 100; // Number of blocks
    IntegralUniform uniform(M, N);
    ofstream WriteData;
    WriteData.open("dataUniform.dat");
        uniform.Run(WriteData);
    WriteData.close();

    IntegralSampling sampling(M, N);
    WriteData.open("dataSampling.dat");
        sampling.Run(WriteData);
    WriteData.close();
    return 0;
}