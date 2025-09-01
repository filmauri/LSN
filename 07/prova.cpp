/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}

double System :: Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){ //Gibbs sampler for Ising
    double d_E = 2.0 * _beta * ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin()) + _H);
    // conditional probability of the spin flip
    (_rnd.Rannyu() < (1./ (1. + exp(-d_E)))) ? _particle(i).setspin(1) : _particle(i).setspin(-1);
    _naccepted++;
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize(bool seed_out){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  
  if(seed_out){
    ifstream Seed("../OUTPUT/seed.out");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);
  }
  else{
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);
  }  

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;                                    
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      // simulation type: 0=MD, 1=MC, 2=Ising MRT^2, 3=Ising Gibbs
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl;
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side;
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "INIT_DELTA" ){
      input >> _init_delta;
      coutf << "INIT_DELTA= " << _init_delta << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();
  if(_sim_type==0) this->initialize_velocities();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: initialize_velocities(){
  double xold, yold, zold;
  if(_restart and (_init_delta == 1)){
    ifstream cinf;
    cinf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> xold >> yold >> zold; // units of coordinates in conf.xyz is _side
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
    // this part was added by me
  } else if(_init_delta == 1){
    // Energy is similar to kT so the velocity is similar to sqrt(kT)
    double fixed_velocity = sqrt(3.0 * _temp); 
    // initialize velocities with a fixed velocity in a random direction
    for (int i = 0; i < _npart; i++) {
      int random_component = int(_rnd.Rannyu(0, 3));
      int direction = int(_rnd.Rannyu(0, 2)) * 2 - 1;
      _particle(i).setvelocity(0, (random_component == 0) ? direction*fixed_velocity : 0.0); 
      _particle(i).setvelocity(1, (random_component == 1) ? direction*fixed_velocity : 0.0); 
      _particle(i).setvelocity(2, (random_component == 2) ? direction*fixed_velocity : 0.0); 
    }

    // set the old positions based on the new velocities
    for (int i = 0; i < _npart; i++) {
      double xold = this->pbc(_particle(i).getposition(0, true) - _particle(i).getvelocity(0) * _delta, 0);
      double yold = this->pbc(_particle(i).getposition(1, true) - _particle(i).getvelocity(1) * _delta, 1);
      double zold = this->pbc(_particle(i).getposition(2, true) - _particle(i).getvelocity(2) * _delta, 2);

      _particle(i).setpositold(0, xold); 
      _particle(i).setpositold(1, yold); 
      _particle(i).setpositold(2, zold); 
    }
  } else if(_restart and _sim_type == 0 ){
    ifstream cinf;
    // initialize velocities from file velocities_.xyz
    cout << "Restarting simulation: Reading velocities from file velocities_.xyz" << endl;
    cinf.open("../OUTPUT/CONFIG/velocities_.xyz");
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & velocities_.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      double vx,vy,vz;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> vx >> vy >> vz; // units of coordinates in conf.xyz is _side
        _particle(i).setvelocity(0, vx);
        _particle(i).setvelocity(1, vy);
        _particle(i).setvelocity(2, vz);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file velocities_.xyz"<< endl;
  } else {
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0);
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // velocity scale factor 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef);
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }

  return;
}

void System :: initialize_properties(string path){ // Initialize data members used for measurement of properties

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  _measure_pofv    = false;

  ifstream input("../INPUT/properties.dat");
  if (input.is_open()){
    std::ostringstream temp_stream;
    temp_stream.precision(2); // 1 decimal place
    temp_stream << std::fixed << _temp;
    string temp_str = temp_stream.str();
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        if(_sim_type > 1){
          cerr << "ERROR: POTENTIAL ENERGY CALCULATION IS NOT APPLICABLE TO THE ISING MODEL SIMULATION." << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutp(path + "potential_energy_" + temp_str + ".dat");
        coutp << "#BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
        _vtail = 8. * M_PI * _rho * (1 / (9. * pow(_r_cut, 9)) - 1 / (3 * pow(_r_cut, 3))); // tail correction for potential energy
      } else if( property == "KINETIC_ENERGY" ){
        if(_sim_type > 1){
          cerr << "ERROR: KINETIC ENERGY CALCULATION IS NOT APPLICABLE TO THE ISING MODEL SIMULATION." << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutk(path + "kinetic_energy_" + temp_str + ".dat");
        coutk << "#BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt(path + "total_energy_" + temp_str + ".dat");
        coutt << "#BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        if(_sim_type > 1){
          cerr << "ERROR: TEMPERATURE CALCULATION IS NOT APPLICABLE TO THE ISING MODEL SIMULATION." << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutte(path + "temperature_" + temp_str + ".dat");
        coutte << "#BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        if(_sim_type > 1){
          cerr << "ERROR: PRESSURE CALCULATION IS NOT APPLICABLE TO THE ISING MODEL SIMULATION." << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpr(path + "pressure_" + temp_str + ".dat");
        coutpr << "#BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 32. * M_PI * _rho * (1 / (9 * pow(_r_cut, 9)) - 1 / (6 * pow(_r_cut, 3))); // tail correction for pressure
      } else if( property == "GOFR" ){
        if(_sim_type > 1){
          cerr << "ERROR: RADIAL DISTRIBUTION FUNCTION CALCULATION IS NOT APPLICABLE TO THE ISING MODEL SIMULATION." << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutgr(path + "gofr_" + temp_str + ".dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << endl;
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "POFV" ){
        if(_sim_type > 0){
          cerr << "ERROR: VELOCITY MODULUS DISTRIBUTION FUNCTION CALCULATION IS NOT APPLICABLE TO THIS SIMULATION TYPE" << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpv(path + "pofv_" + temp_str + ".dat");
        coutpv << "#VELOCITY:     AVE_POFV:     PROGR_AVE_POVF:        ERROR:" << endl;
        coutpv.close();
        input>>_n_bins_v;
        _nprop += _n_bins_v;
        _bin_size_v = 4.0*sqrt(_temp)/(double)_n_bins_v; 
        _measure_pofv = true;
        _index_pofv = index_property;
        index_property += _n_bins_v;
      } else if( property == "MAGNETIZATION" ){
        if(_sim_type < 2){
          cerr << "ERROR: MAGNETIZATION CALCULATION IS NOT APPLICABLE TO LJ SIMULATION" << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpr(path + "magnetization_" + temp_str + ".dat");
        coutpr << "#BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr(path + "specific_heat_" + temp_str + ".dat");
        coutpr << "#BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        if(_sim_type < 2){
          cerr << "ERROR: SUSCEPTIBILITY CALCULATION IS NOT APPLICABLE TO LJ SIMULATION" << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpr(path + "susceptibility_" + temp_str + ".dat");
        coutpr << "#BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
        // when we encounter ENDPROPERTIES the parsing stops 
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System:: reset_properties(){ // Reset properties to zero
  _measurement.zeros();
  _average.zeros();
  _block_av.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System :: finalize(string path){
  if (_measure_gofr) {
    ofstream coutf(path + "gofr_final.dat");
    coutf << "#DISTANCE:     AVE_GOFR:        ERROR:" << endl;
    for (int i = 0; i < _n_bins; i++) {
      double r = (i + 0.5) * _bin_size; 
      double sum_average = _global_av(_index_gofr + i);
      double sum_ave2 = _global_av2(_index_gofr + i);

      coutf << setw(12) << r
            << setw(12) << sum_average / (double(_nblocks)) 
            << setw(12) << this->error(sum_average, sum_ave2, _nblocks) << endl;
    }
    coutf.close();
  }
  // write the final configuration and velocities
  this->write_configuration();
  this->write_velocities();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  "
              << setprecision(17) << _particle(i).getposition(0,false)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,false)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,false)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
    coutf.close();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// write velocity configuration as .xyz file in directory ../OUTPUT/CONFIG/

void System :: write_velocities(){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/velocities_.xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getvelocity(0)          // vx
            << setw(16) << _particle(i).getvelocity(1)          // vy
            << setw(16) << _particle(i).getvelocity(2) << endl; // vz
    }
  } else cerr << "PROBLEM: Unable to open velocity.xyz" << endl;
  coutf.close();
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}


// Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System :: read_configuration(){
  ifstream cinf;
  // read input fcc
  if(_sim_type <= 1 and !_restart) {
  cinf.open("../INPUT/CONFIG/config.fcc");
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cout << "ncoord = " << ncoord << " _npart = " << _npart << endl;
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.fcc"<< endl;
  }
  
  // restart simulation type 0 or 1
  if(_sim_type <= 1 and _restart){
    cinf.open("../OUTPUT/CONFIG/config.xyz");
    cout << "Restarting simulation: Reading configuration from file config.xyz" << endl;
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      double x,y,z;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
        _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
        _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
        _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
    cinf.close();
    cinf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & conf-1.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      double xold,yold,zold;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> xold >> yold >> zold; // units of coordinates in conf.xyz is _side
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
  }
  // read input spin
  if(_sim_type > 1 and !_restart){
    cinf.open("../INPUT/CONFIG/config.spin");
    if(cinf.is_open()){
      int spin;
      for(int i=0; i<_npart; i++){
        cinf >> spin;
        _particle(i).setspin(spin);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file config.spin"<< endl;
    cinf.close();
  }
  // restart simulation type 2 or 3
  if(_restart and _sim_type > 1){
    if(cinf.is_open()){
      cinf.open("../OUTPUT/CONFIG/config.spin");
      int spin;
      for(int i=0; i<_npart; i++){
        cinf >> spin;
        _particle(i).setspin(spin);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file config.spin"<< endl;
    cinf.close();
  } 
  return;
}

void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}


void System :: measure(){ // Measure properties
  _measurement.zeros();
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;
  double magnetization=0.0;
  double virial=0.0;
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        // GOFR ... TO BE FIXED IN EXERCISE 7
        if(_measure_gofr){
          bin = int(dr/_bin_size); //bin index corresponding to distance r
          if(bin < _n_bins) _measurement(_index_gofr + bin) += 2.0; // factor 2 for symmetry
        }
        if(dr < _r_cut){
          if(_measure_penergy)  penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          if(_measure_pressure) virial       += 1.0/pow(dr,12) - 0.5/pow(dr,6); // PRESSURE
        }
      }
    }
  }
  // POFV ... TO BE FIXED IN EXERCISE 4, modified 
  if(_measure_pofv){
    for(int i=0; i<_npart; i++){
      double v = sqrt( _particle(i).getvelocity(0) * _particle(i).getvelocity(0) +
                       _particle(i).getvelocity(1) * _particle(i).getvelocity(1) +
                       _particle(i).getvelocity(2) * _particle(i).getvelocity(2) );
      bin = int(v/_bin_size_v); //bin index corresponding to velocity v
      if(bin < _n_bins_v) _measurement(_index_pofv + bin) += 1.0;
    }
  }
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if ((_measure_penergy or _measure_tenergy) and _sim_type < 2){ // if i want to evaluate total energy for LJ i need penergy
  //if (_measure_penergy){ come il prof
  //if (_measure_penergy or _measure_tenergy or _measure_cv){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    if (_measure_penergy) _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  //if ((_measure_kenergy or _measure_tenergy or _measure_temp or _measure_pressure) and _sim_type < 2){
  if (_measure_kenergy or _measure_temp){ //come il prof
  //if (_measure_kenergy or _measure_tenergy or _measure_temp or _measure_pressure or _measure_cv){
    if(_sim_type == 0){
      for (int i=0; i<_npart; i++){
        kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() );
      } 
      kenergy_temp /= double(_npart);
    }
    if(_sim_type == 1) kenergy_temp += 0.5 * 3. * _temp; 
    
    if (_measure_kenergy) _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy or _measure_cv){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
    }
    if(_measure_tenergy) _measurement(_index_tenergy) = tenergy_temp;
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){ // to compute temperature i need kenergy, NVE will give T, NVT 0
  //if (_measure_temp or _measure_kenergy){
    _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    if(_sim_type == 0){
      _measurement[_index_pressure] = _rho * (2.0/3.0) * kenergy_temp + (_ptail*_npart + 48.0*virial/3.0)/(_volume);
    }
    if(_sim_type == 1){
      _measurement[_index_pressure] = _rho * _temp + (_ptail*_npart + 48.0*virial/3.0)/(_volume);
    }
  } 
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if(_measure_magnet or _measure_chi){
    for(int i=0; i<_npart; i++){
      magnetization += double(_particle(i).getspin());
    }
    if (_measure_magnet) _measurement(_index_magnet) = magnetization;
    
  }
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if (_measure_cv){
    _measurement(_index_cv) = pow(tenergy_temp * double(_npart),2);
  }
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  if(_measure_chi ){
    //double chi_temp = _measurement(_index_magnet)*_measurement(_index_magnet)*_beta;
    double chi_temp = magnetization*magnetization * _beta / double(_npart);
    _measurement(_index_chi) = chi_temp;
  }

  _block_av += _measurement; //Update block accumulators
 
  return;
}

void System :: averages(int blk, string path){

  ofstream coutf;
  ostringstream temp_stream;
  temp_stream.precision(2); 
  temp_stream << fixed << _temp;

  double average, sum_average, sum_ave2;

  if(_measure_cv){ 
    // subtract <H>^2 to the block average and multiply by beta^2
    _block_av(_index_cv) -= pow(_block_av(_index_tenergy) * double(_npart), 2) / double(_nsteps);
    _block_av(_index_cv) *= _beta * _beta;
  }

  if (_measure_gofr) {
    for (int i = 0; i < _n_bins; i++) {
      //double r = (i + 0.5) * _bin_size; 
      double shell_volume = (4.0 / 3.0) * M_PI * (pow(static_cast<double>(i + 1) * _bin_size, 3) - pow(static_cast<double>(i) * _bin_size, 3)); 
      double norm = shell_volume * _rho * double(_npart); 
      _block_av(_index_gofr + i) /= norm; 
    }
  }

  _average = _block_av / double(_nsteps);

  if(_measure_magnet){
    _average(_index_magnet) /= double(_npart);
  }

  if(_measure_cv){
    _average(_index_cv) /= double(_npart);
  }

  _global_av  += _average;
  _global_av2 += _average % _average;

  if (_measure_penergy){
    coutf.open(path + "/potential_energy_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if (_measure_kenergy){
    coutf.open(path + "/kinetic_energy_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if (_measure_tenergy){
    coutf.open(path + "/total_energy_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if (_measure_temp){
    coutf.open(path + "/temperature_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if (_measure_pressure){
    coutf.open(path + "/pressure_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if(_measure_pofv){
    coutf.open(path + "/pofv_" + temp_stream.str() + ".dat", ios::app);
    for(int i=0; i<_n_bins_v; i++){
      average  = _average(_index_pofv + i);
      sum_average = _global_av(_index_pofv + i);
      sum_ave2 = _global_av2(_index_pofv + i);
      coutf << setw(12) << (i+0.5)*_bin_size_v
            << setw(12) << average/double(_npart)
            << setw(12) << sum_average/(double(blk)*double(_npart))
            << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    }
    coutf.close();
  }

  if(_measure_magnet){
    coutf.open(path + "/magnetization_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if(_measure_cv){
    coutf.open(path + "/specific_heat_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_cv);
    sum_average = _global_av(_index_cv);
    sum_ave2 = _global_av2(_index_cv);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if(_measure_chi){
    coutf.open(path + "/susceptibility_" + temp_stream.str() + ".dat", ios::app);
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }

  if (_measure_gofr) {
    ofstream coutf(path + "/gofr_" + temp_stream.str() + ".dat", ios::app);
    for (int i = 0; i < _n_bins; i++) {
        double r = (i + 0.5) * _bin_size;
        double average = _average(_index_gofr + i);
        double sum_average = _global_av(_index_gofr + i);
        double sum_ave2 = _global_av2(_index_gofr + i);

        coutf << setw(12) << r
              << setw(12) << sum_average / double(blk) 
              << setw(12) << this->error(sum_average, sum_ave2, blk)  << endl; 
    }
    coutf.close();
  }

  
  double fraction;

  if(_sim_type > 0){
    coutf.open(path + "/acceptance_" + temp_stream.str() + ".dat", ios::app);
    if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
    else fraction = 0.0; 
    coutf << setw(12) << blk << setw(12) << fraction << endl;
    coutf.close();
  }
  return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
}

void System :: invert_velocities(){
  for(int i=0; i<_npart; i++){
    _particle(i).setvelocity(0, -_particle(i).getvelocity(0));
    _particle(i).setvelocity(1, -_particle(i).getvelocity(1));
    _particle(i).setvelocity(2, -_particle(i).getvelocity(2));
  }
  return;
}

void System :: set_temp(double temp){
  _temp = temp;
  _beta = 1.0/_temp;
  return;
}

int System :: get_sim_type(){
  return _sim_type;
}

int System :: get_naccepted(){
  return _naccepted;
}

int System :: get_nattempts(){
  return _nattempts;
}

bool System :: get_restart(){
  return _restart;
}