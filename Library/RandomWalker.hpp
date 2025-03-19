#ifndef __Random_Walker__
#define __Random_Walker__

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "PRNG/random.hpp"
#include "library.hpp"
using namespace std; 

class Walker{
    public:
        //constructor
        Walker(){
            position = vector<double>(3, 0.0);
            oldPosition = vector<double>(3, 0.0);
        }
        ~Walker() {;}

        void SetPosition(double x_, double y_, double z_){
            x = x_;
            y = y_;
            z = z_;
        }
        vector <double> getPosition(){
            return position;
        }

        void latticeStep(int y){
            //int y = rnd.DiscreteMovement();
            //cout << "numero casuale "<< y << endl;
            position[y/2] += ((y % 2 == 0) ? 1. : -1.);
        }
        void continuousStep(double theta, double phi){
            position[0] += sin(theta)*cos(phi);
            position[1] += sin(theta)*sin(phi);
            position[2] += cos(theta);
        }

        double getDistance(){
            return (position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
        }
        double getX(){return position[0];}
        double getY(){return position[1];}
        double getZ(){return position[2];}
        private:
            double x, y, z;
            //Random rnd;
            vector <double> position;
            vector <double> oldPosition;
};

#endif