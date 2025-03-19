#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "../../Library/PRNG/random.hpp"
#include "../../Library/blockAverage/blockAverage.hpp"
#include "../../Library/library.hpp"
#include "../../Library/RandomWalker.hpp"

using namespace std;
int main()
{   
    double M = 10000; //number of throws
    double N = 100; //number of blocks
    double L = M/N; //number of throws in each block
    double length = 100; //lenght of the random walk

    ofstream WriteDiscrete;
    WriteDiscrete.open("discrete.dat");
    if (!WriteDiscrete.is_open())
        cerr << "PROBLEM: Unable to open discrete.dat" << endl;

    ofstream WriteDiscreteWalk;
    WriteDiscreteWalk.open("discreteWalk.dat");
    if (!WriteDiscreteWalk.is_open())
        cerr << "PROBLEM: Unable to open discreteWalk.dat" << endl;

    ofstream WriteContinuous;
    WriteContinuous.open("continuous.dat");
    if (!WriteContinuous.is_open())
        cerr << "PROBLEM: Unable to open continuous.dat" << endl;
    
    ofstream WriteContinuousWalk;
    WriteContinuousWalk.open("continuousWalk.dat");
    if (!WriteContinuousWalk.is_open())
        cerr << "PROBLEM: Unable to open continuousWalk.dat" << endl;

    Random rnd("../../Library/PRNG/");
    Walker walker;
    vector <Walker> discrete(M, walker);
    vector <Walker> continuous(M, walker);

    double blockDiscrete = 0, blockDiscrete2 = 0;
    double blockDiscreteProg = 0, blockDiscreteProg2 = 0;
    double errorDiscrete = 0;

    double blockContinuous = 0, blockContinuous2 = 0;
    double blockContinuousProg = 0, blockContinuousProg2 = 0;
    double errorContinuous = 0;

    for (int l = 0; l < length; l++){       //ad ogni iterazione viene simulato un passo per tutti i rw
        for (int b = 0; b < N; b++){        //ciclo sui blocchi
            for (int i = 0; i < L; i++){    //ciclo sulle simulazioni
                discrete[L * b + i].latticeStep(rnd.DiscreteMovement());
                blockDiscrete2 += discrete[L * b + i].getDistance();

                continuous[L * b + i].continuousStep(rnd.Angle(), rnd.Angle()/2);
                blockContinuous2 += continuous[L * b + i].getDistance();
            }
            blockDiscrete2 /= L; //media della distanza quadra sul blocco
            blockDiscrete = sqrt(blockDiscrete2); //radice della media della distanza quadra sul blocco
            blockDiscreteProg += blockDiscrete;
            blockDiscreteProg2 += blockDiscrete2;
            blockDiscrete = 0, blockDiscrete2 = 0;

            blockContinuous2 /= L;
            blockContinuous = sqrt(blockContinuous2);
            blockContinuousProg += blockContinuous;
            blockContinuousProg2 += blockContinuous2;
            blockContinuous = 0, blockContinuous2 = 0;

        }
        blockDiscreteProg /= N; //media progressiva
        blockDiscreteProg2 /= N;
        errorDiscrete = computeError(blockDiscreteProg, blockDiscreteProg2, N);
        WriteDiscrete << l << " " << blockDiscreteProg << " " << errorDiscrete << endl;
        WriteDiscreteWalk << l << " " << discrete[0].getX() << " " <<discrete[0].getY() << " " << discrete[0].getZ() << endl;
        blockDiscreteProg = 0, blockDiscreteProg2 = 0;
        errorDiscrete = 0;

        blockContinuousProg /= N;
        blockContinuousProg2 /= N;
        errorContinuous = computeError(blockContinuousProg, blockContinuousProg2, N);
        WriteContinuous << l << " " << blockContinuousProg << " " << errorContinuous << endl;
        WriteContinuousWalk << l << " " << continuous[0].getX() << " " <<continuous[0].getY() << " " << continuous[0].getZ() << endl;
        blockContinuousProg = 0, blockContinuousProg2 = 0;
        errorContinuous = 0;
    }
    
    return 0;
}