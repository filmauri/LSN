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
    double M = 1000; //number of throws
    double N = 100; //number of blocks
    double L = M/N; //number of throws in each block
    double length = 100000; //lenght of the random walk

    ofstream WriteContinuous;
    WriteContinuous.open("continuous.dat");
    if (!WriteContinuous.is_open())
        cerr << "PROBLEM: Unable to open continuous.dat" << endl;

    ofstream WriteContinuousWrong;
    WriteContinuousWrong.open("continuousWrong.dat");
    if (!WriteContinuousWrong.is_open())
        cerr << "PROBLEM: Unable to open continuousWrong.dat" << endl;
    
    ofstream WriteContinuousWalk;
    WriteContinuousWalk.open("continuousWalk.dat");
    if (!WriteContinuousWalk.is_open())
        cerr << "PROBLEM: Unable to open continuousWalk.dat" << endl;

    ofstream WriteContinuousWrongWalk;
    WriteContinuousWrongWalk.open("continuousWrongWalk.dat");
    if (!WriteContinuousWrongWalk.is_open())
        cerr << "PROBLEM: Unable to open continuousWrongWalk.dat" << endl;

    Random rnd("../../Library/PRNG/");
    Walker walker;
    vector <Walker> continuous(M, walker);
    vector <Walker> wrong_continuous(M, walker);

    double blockContinuous = 0, blockContinuous2 = 0;
    double blockContinuousProg = 0, blockContinuousProg2 = 0;
    double errorContinuous = 0;

    double blockWrongContinuous = 0, blockWrongContinuous2 = 0;
    double blockWrongContinuousProg = 0, blockWrongContinuousProg2 = 0;
    double errorWrongContinuous = 0;


    for (int l = 0; l < length; l++){       //ad ogni iterazione viene simulato un passo per tutti i rw
        for (int b = 0; b < N; b++){        //ciclo sui blocchi
            for (int i = 0; i < L; i++){    //ciclo sulle simulazioni

                continuous[L * b + i].continuousStep(rnd.Angle(), rnd.phi_spherical(rnd.Rannyu()));
                blockContinuous2 += continuous[L * b + i].getDistanceXY();

                wrong_continuous[L * b + i].continuousStep(rnd.Angle(), rnd.Angle()/2);
                blockWrongContinuous2 += wrong_continuous[L * b + i].getDistanceXY();
            }

            blockContinuous2 /= L;
            blockContinuous = sqrt(blockContinuous2);
            blockContinuousProg += blockContinuous;
            blockContinuousProg2 += blockContinuous2;
            blockContinuous = 0, blockContinuous2 = 0;

            blockWrongContinuous2 /= L;
            blockWrongContinuous = sqrt(blockWrongContinuous2);
            blockWrongContinuousProg += blockWrongContinuous;
            blockWrongContinuousProg2 += blockWrongContinuous2;
            blockWrongContinuous = 0, blockWrongContinuous2 = 0;

        }
        blockContinuousProg /= N;
        blockContinuousProg2 /= N;
        errorContinuous = computeError(blockContinuousProg, blockContinuousProg2, N);
        WriteContinuous << l << " " << blockContinuousProg << " " << errorContinuous << endl;
        WriteContinuousWalk << l << " " << continuous[0].getX() << " " <<continuous[0].getY() << " " << continuous[0].getZ() << endl;
        blockContinuousProg = 0, blockContinuousProg2 = 0;
        errorContinuous = 0;

        blockWrongContinuousProg /= N;
        blockWrongContinuousProg2 /= N;
        errorWrongContinuous = computeError(blockWrongContinuousProg, blockWrongContinuousProg2, N);
        WriteContinuousWrongWalk << l << " " << wrong_continuous[0].getX() << " " <<wrong_continuous[0].getY() << " " << wrong_continuous[0].getZ() << endl;
        WriteContinuousWrong << l << " " << blockWrongContinuousProg << " " << errorWrongContinuous << endl;
        blockWrongContinuousProg = 0, blockWrongContinuousProg2 = 0;
        errorWrongContinuous = 0;
    }
    
    return 0;
}