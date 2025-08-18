#ifndef __CITIES__
#define __CITIES__
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;

class City {
public:
    City(){
        _coordinates.resize(_dim);
    };
    void setLocation(arma::vec coordinates) 
    {
        if (_coordinates.size() != _dim) {
            cerr << "Error: Dimension mismatch in setLocation." << endl;
            return;
        }
        _coordinates = coordinates;
    }
    
    arma::vec getCoordinates(){return _coordinates;}
    double getX(){return _coordinates(0);}
    double getY(){return _coordinates(1);}
    
private:
    arma::vec _coordinates; // Vector to store coordinates of cities
    const int _dim = 2; // 2D space 
};


#endif