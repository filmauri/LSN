#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "../../Library/PRNG/random.hpp"
#include "../../Library/blockAverage/blockAverage.hpp"
#include "../../Library/library.hpp"

using namespace std;
int main(){
    Random rnd("../../Library/PRNG/");
    int M = 10000; // Number of throws
    int N = 100; // Number of blocks

    callOption call(M, N);
    ofstream WriteData;
    WriteData.open("dataCall1step.dat");
        call.Run(WriteData);
    WriteData.close();

    callOptionStep callStep(M, N);
    ofstream WriteDataCallStep;
    WriteDataCallStep.open("dataCallSteps.dat");
        callStep.Run(WriteDataCallStep);
    WriteDataCallStep.close();

    putOption put(M, N);
    ofstream WriteDataPut;
    WriteDataPut.open("dataPut1step.dat");
        put.Run(WriteDataPut);
    WriteDataPut.close();

    putOptionStep putStep(M, N);
    ofstream WriteDataPutStep;
    WriteDataPutStep.open("dataPutSteps.dat");
        putStep.Run(WriteDataPutStep);
    WriteDataPutStep.close();





    return 0;
}