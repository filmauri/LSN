#ifndef __GENETICLIBRARY__
#define __GENETICLIBRARY__
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>
#include "../../../Library/PRNG/random.hpp"
#include "../../../Library/library.hpp"
#include "cities.hpp"
#include "mpi.h"
using namespace std;
using namespace arma;

void initializeRandomDirectly(Random& rnd, int rank) {
    // Leggi Primes
    ifstream Primes("../../../Library/PRNG/Primes");
    int p1, p2;
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
        Primes.close();
    }
    
    // Leggi seed base
    ifstream input("../../../Library/PRNG/seed.in");
    string property;
    int baseSeed[4];
    
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> baseSeed[0] >> baseSeed[1] >> baseSeed[2] >> baseSeed[3];
                break;
            }
        }
        input.close();
    }
    
    // Modifica i semi in base al rank
    int seeds[4];
    seeds[0] = (baseSeed[0] + rank * 1237) % 4096;
    seeds[1] = (baseSeed[1] + rank * 2347) % 4096; 
    seeds[2] = (baseSeed[2] + rank * 3457) % 4096;
    seeds[3] = (baseSeed[3] + rank * 4567) % 4096;
    
    // Usa SetRandom per bypassare il file
    rnd.SetRandom(seeds, p1, p2);
    
    cout << "Rank " << rank << ": Direct seeds [" << seeds[0] << ", " 
         << seeds[1] << ", " << seeds[2] << ", " << seeds[3] << "]" << endl;
}

class TSP{
    
public: 
    //reading the input file and initializing the parameters
    void initialize(int rank){
        _rnd = Random("../../../Library/PRNG/");
        initializeRandomDirectly(_rnd, rank);
        cout << "Rank " << rank << ": Random number generator initialized." << _rnd.Rannyu() <<endl;
        ifstream input("../INPUT/input.dat");
        _rank = rank;
        if (!input.is_open()) {
            cerr << "Error opening input file." << endl;
            exit(EXIT_FAILURE);
        }
        ofstream coutf("../OUTPUT/output.txt");
        string property;
        while(!input.eof()){
            input >> property;
            if (property == "TYPE") {
                input >> _type;
                if (_type > 2) {
                    cerr << "Error: Unsupported problem type.\n 0 for Cities on a circumference, 1 for cities scattered in a square" << endl;
                    exit(EXIT_FAILURE);
                }
                if (_type == 0) {
                    coutf.open("output.txt");
                    coutf << "Problem type: Cities on a Circle" << endl;
                } else if (_type == 1) {
                    coutf.open("output.txt");
                    coutf << "Problem type: Cities scattered in a square" << endl;
                } else if(_type == 2) {
                    coutf.open("output.txt");
                    coutf << "Problem type: Il giro d'Italia" << endl;
                }
            } else if(property == "MIGRATION_STEP"){
                input >> _migrationStep;
                coutf.open("output.txt");
                coutf << "Migration step: " << _migrationStep << endl;

            }else if (property == "NORM_ORDER"){
                input >> _normOrder;
                if (_normOrder >= 1) 
                    coutf << "Norm order" << _normOrder << endl;
                else{
                    cerr << "Error: Norm order must be at least 1." << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (property == "POWER") {
                input >> _power;
                if (_power >= 1) {
                    coutf << "Power: " << _power << endl;
                }
                else {
                    cerr << "Error: Power must be at least 1." << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (property == "N_CITIES") {
                input >> _nCities;
                if (_nCities > 2) {
                    coutf << "Number of cities" << _nCities << endl;
                }
                else{
                    cerr << "Error: Number of cities must be at least 2." << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (property == "N_INDIVIDUALS")
            {
                input >> _nIndividuals;
                if (_nIndividuals > 0) {
                    coutf << "Number of individuals in the population: " << _nIndividuals << endl;
                }
                else{
                    cerr << "Error: Number of individuals must be at least 1." << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (property == "N_GENERATIONS") {
                input >> _nGenerations;
                if (_nGenerations > 0) {
                    coutf << "Number of generations: " << _nGenerations << endl;
                }
                else{
                    cerr << "Error: Number of generations must be at least 1." << endl;
                    exit(EXIT_FAILURE);
                }
            } else if(property == "PROB_MUTATIONS"){
                int size;
                input >> size;
                _probMutations.set_size(size);
                for (int i = 0; i < size; i++) {
                    input >> _probMutations(i);
                    if (_probMutations(i) < 0 || _probMutations(i) > 1) {
                        cerr << "Error: Mutation probability must be between 0 and 1." << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                coutf << "Mutation probabilities initialized." << endl;
                
                
            }else if (property == "ENDINPUT") {
                _loss.set_size(_nIndividuals);
                coutf << "System initialized!" << std::endl;
                coutf.close();
                break;
            } else {
                cerr << "Error: Unknown property in input file." << endl;
                exit(EXIT_FAILURE);
            }
            
        }
        _cities.set_size(_nCities);
        if (_rank == 0) {
        cout << "Rank " << _rank << ": Initializing cities..." << endl;
        this->initialize_cities();
        cout << "Rank " << _rank << ": Generating TSP population..." << endl;
        this->generate_tsp_population();
        cout << "Rank " << _rank << ": Checking generated population" << endl;
        this->checkPopulation();
        this->checkStartingPos();
        cout << "Rank " << _rank << ": Population check completed." << endl;
    } 
    }
    
    void initialize_cities(){
       // _cities.set_size(_nCities);
        arma::vec pos;
        pos.set_size(_dimension);
        ifstream citiesFile("../INPUT/cap_prov_ita.dat");
        if (!citiesFile.is_open()) {
            cerr << "Error opening cities file." << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < _nCities; i++){
            if (_type == 0) { // Circle
                double angle = _rnd.Angle();
                pos(0) = cos(angle);
                pos(1) = sin(angle);
                _cities(i).setLocation(pos);
            } else if (_type == 1) { // Scatter in a square
                pos(0) = _rnd.Rannyu(-1, 1);
                pos(1) = _rnd.Rannyu(-1, 1);
                _cities(i).setLocation(pos);
            }
            else if (_type == 2) { // Giro d'Italia
                double x, y;
                citiesFile >> x >> y; // Read first city
                pos(0) = x;
                pos(1) = y;
                _cities(i).setLocation(pos);
                if (i == 0) {
                    cout << "Giro d'Italia cities initialized." << endl;
                }
            } else {
                cerr << "Error: Unsupported problem type." << endl;
                exit(EXIT_FAILURE);
            } 
        }
        ofstream coutf("../OUTPUT/output.txt", ios::app);
        for (int i = 0; i < _nCities; i++) {
            coutf << "City " << i << ": (" << _cities(i).getX() << ", " << _cities(i).getY() << ")" << endl;
        }
        coutf << "Cities initialized." << endl;
        coutf.close();
    }

void generate_tsp_population() {
    // PROBLEMA: Usa parametri invece di variabili membro
    _population.set_size(_nCities, _nIndividuals);
    arma::vec base_remaining = arma::regspace(1, _nCities - 1);
    
    for (int i = 0; i < _nIndividuals; ++i) {
        _population(0, i) = 0;
        arma::vec shuffled = arma::shuffle(base_remaining);
        _population.submat(1, i, _nCities - 1, i) = shuffled;
    }
}

void checkPopulation() {
    if (_population.n_rows != _nCities || _population.n_cols != _nIndividuals) {
        cerr << "Error: Population size mismatch." << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < _nIndividuals; ++i) {
        arma::vec individual = _population.col(i);
        for (int j = 0; j < _nCities; ++j) {
            if (arma::sum(individual == j) != 1) {
                cerr << "Error: Individual " << i << " does not contain all cities." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}

void checkStartingPos(){
    int invalid_count = 0;
    for (int i = 0; i < _population.n_cols; ++i) 
        if (_population(0, i) != 0) 
            invalid_count++;
        
    if (invalid_count != 0) 
        exit(EXIT_FAILURE);
}

void printPopulation(){
    cout << "Population of TSP: " << _nIndividuals << endl;
    
    for (int i = 0; i < _nIndividuals; i++) {
        cout << "Individual " << i + 1 << ": ";
        for (int j = 0; j < _nCities; j++) {
            cout << _population(j, i) << " ";
        }
        cout << endl;
    }
}

double distance(int city, int individualIndex) {
    int city1_index = _population(city, individualIndex);
    int city2_index = _population(boundary_condition(city + 1), individualIndex);
    
    arma::vec city1_pos = _cities(city1_index).getCoordinates();
    arma::vec city2_pos = _cities(city2_index).getCoordinates();
    
    double dx = city1_pos(0) - city2_pos(0);
    double dy = city1_pos(1) - city2_pos(1);
    
    if (_normOrder == 1) {
        return sqrt(dx*dx + dy*dy);  // Distanza euclidea
    } else if (_normOrder == 2)     
        return dx*dx + dy*dy; // Distanza quadratica
    else 
        return 1;
}

int boundary_condition(int i_city) {
    if (i_city < 0) {
        return _nCities - 1; // Torna all'ultima città se si esce a sinistra
    } else if (i_city >= _nCities) {
        return 0; // Torna alla prima città se si esce a destra
    }
    return i_city; // Ritorna l'indice della città se è valido
}

double loss_travel(int individualIndex) {
    double loss = 0.0;
    for (int cityIndex = 0; cityIndex < _nCities; cityIndex++) {
        loss += distance(cityIndex, individualIndex);
    }
    return loss;
}

void loss(){
    for (int travelIndex = 0; travelIndex < _nIndividuals; travelIndex++) {
        double individualLoss = loss_travel(travelIndex);
        if (individualLoss < 0) {
            cerr << "Error: Negative loss encountered for individual " << travelIndex << "." << endl;
            exit(EXIT_FAILURE);
        }
        _loss(travelIndex) = individualLoss;
    }   
    this->sortPopulationbyLoss();
}

void sortPopulationbyLoss() { 
    arma::uvec sortedIndices = arma::sort_index(_loss);
    _population = _population.cols(sortedIndices);
    _loss = _loss(sortedIndices);
}

void bestTravel(int gen) {
    ofstream coutf("../OUTPUT/output.txt", ios::app);
    if (gen % 50 == 0) {
        coutf << "Best travel in generation " << gen << ": ";
        for (int i = 0; i < _nCities; i++) {
            coutf << _population(i, 0) << " ";
        }
        coutf << "Loss: " << _loss(0) << endl;
        coutf << endl;
    }
    coutf.close();
}
int selection() {
    int index;
    index = static_cast<int>(_nIndividuals * (pow(_rnd.Rannyu(), _power)));
    return index;
}

void swap (int a, int b){
    int t = a;
	a = b;
	b = t;
}

void pairPermutationWithFixedStart(int individual) {
    
    int city1 = 1 + int(_rnd.Rannyu() * (_nCities - 1));
    int city2;
    do {
        city2 = 1 + int(_rnd.Rannyu() * (_nCities - 1));
    } while (city1 == city2);

    double temp = _population(city1, individual);
    _population(city1, individual) = _population(city2, individual);
    _population(city2, individual) = temp;
    
    this->checkPopulation();
    this->checkStartingPos();
}

void shift(int travelIndex) {
    int shift = int(_rnd.Rannyu() * (_nCities - 2)) + 1; // Random shift between 1 and n_cities-2
    
    // Explicit conversion from double to int - correct syntax
    arma::Col<int> temp_travel = arma::conv_to<arma::Col<int>>::from(_population.col(travelIndex));
    
    arma::Col<int> route_without_start = temp_travel.subvec(1, _nCities - 1);
    
    // Perform circular shift on the remaining cities
    for (int i = 0; i < shift; i++) {
        int last_city = route_without_start(route_without_start.n_elem - 1); // Get last city
        route_without_start.shed_row(route_without_start.n_elem - 1);         // Remove last city
        route_without_start.insert_rows(0, arma::Col<int>({last_city}));      // Insert at beginning
    }
    
    // Reconstruct the complete travel with the first city still at position 0
    temp_travel.subvec(1, _nCities - 1) = route_without_start;
    
    // Convert back to double if needed and assign
    _population.col(travelIndex) = arma::conv_to<arma::Col<double>>::from(temp_travel);

    this->checkPopulation();
    this->checkStartingPos();
}

void blockPermutation(int travelIndex) {
    int n_index = int(_rnd.Rannyu(1, _nCities/4)); // Block size between 1 and nCities/4
    int block1_start = int(_rnd.Rannyu(1, _nCities/2)); // Start index for block 1
    int block2_start = _nCities/2; 
    do {
        block2_start = int(_rnd.Rannyu(1, _nCities - n_index)); // Start index for block 2
    }while(block2_start - block1_start <= n_index);

    arma::Col<int> temp_travel = arma::conv_to<arma::Col<int>>::from(_population.col(travelIndex));
    
    // Estrai i due blocchi
    arma::Col<int> block1 = temp_travel.subvec(block1_start, block1_start + n_index);
    arma::Col<int> block2 = temp_travel.subvec(block2_start, block2_start + n_index);
    
    // Scambia i blocchi
    temp_travel.subvec(block1_start, block1_start + n_index) = block2;
    temp_travel.subvec(block2_start, block2_start + n_index) = block1;
    
    // Riconverti e aggiorna la popolazione
    _population.col(travelIndex) = arma::conv_to<arma::colvec>::from(temp_travel);
    
    // Verifica la validità della soluzione
    this->checkPopulation();
    this->checkStartingPos();
}

void Inversion(int travelIndex) {
    // Select two random indices to define the block
    int start = int(_rnd.Rannyu() * (_nCities - 2)) + 1;  // Range: 1 to nCities-2
    int end = int(_rnd.Rannyu() * (_nCities - start)) + start;  // Range: start to nCities-1
    if (start > end) {
        swap(start, end);
    }

    arma::Col<int> temp_travel = arma::conv_to<arma::Col<int>>::from(_population.col(travelIndex));
    arma::Col<int> block = temp_travel.subvec(start, end);
    block = arma::flipud(block);
    
    temp_travel.subvec(start, end) = block;
    _population.col(travelIndex) = arma::conv_to<arma::colvec>::from(temp_travel);
    
    this->checkPopulation();
    this->checkStartingPos();
}

void mutation(int travelIndex) {
    if (_rnd.Rannyu() < _probMutations(0)) {
        pairPermutationWithFixedStart(travelIndex);
    } else if (_rnd.Rannyu() < _probMutations(1)) {
        shift(travelIndex);
    } else if (_rnd.Rannyu() < _probMutations(2)) {
        blockPermutation(travelIndex);
    } else if (_rnd.Rannyu() < _probMutations(3)) {
        Inversion(travelIndex);
    } else {
        // No mutation
        return;
    }
} 

void crossOver(int travelIndex, int index1, int index2) {
    int cutPoint = 1 + _rnd.Rannyu(0, _nCities - 2); // Taglio tra posizione 1 e _nCities-2
    
    arma::Col<int> parent1 = arma::conv_to<arma::Col<int>>::from(_population.col(index1));
    arma::Col<int> parent2 = arma::conv_to<arma::Col<int>>::from(_population.col(index2));

    // 3. Crea i due offspring
    arma::Col<int> offspring1(_nCities);
    arma::Col<int> offspring2(_nCities);
    
    // 4. 5. Conserva la prima parte di entrambi i genitori e crea liste delle città utilizzate
    std::vector<bool> used1(_nCities, false);
    std::vector<bool> used2(_nCities, false);
    for (int i = 0; i < cutPoint; i++) {
        offspring1(i) = parent1(i);
        offspring2(i) = parent2(i);
        used1[offspring1(i)] = true;
        used2[offspring2(i)] = true;
    }

    // 6. Completa offspring1 con le città mancanti nell'ordine di parent2
    int pos1 = cutPoint;
    for (int i = 0; i < _nCities && pos1 < _nCities; i++) {
        int city = parent2(i);
        if (!used1[city]) {
            offspring1(pos1) = city;
            used1[city] = true;
            pos1++;
        }
    }
    
    _population.col(travelIndex) = arma::conv_to<arma::colvec>::from(offspring1);

    // 9. Verifica la validità
    this->checkPopulation();
    this->checkStartingPos();
}

//perform one generation of the genetic algorithm
void evolution(int num) {
    this->loss(); // Ensure the population is sorted by loss
    
    for (int i = 0; i < _nIndividuals; i++) {
        int parent1 = this->selection();
        int parent2; 
        do {
            parent2 = this->selection();
        }while(parent2 == parent1); 

        if (_rnd.Rannyu() < _probMutations(4)) { 
        this->crossOver(i, parent1, parent2);
        }

        this->mutation(i);    
        this->loss();
    }
    
    this->checkPopulation();
    this->checkStartingPos();

    if (num % 10 == 0) {
        this->partialFinalization(num);
    }
} 


int getNGenerations() const {
    return _nGenerations;
}

void finalize() {
    ofstream coutf("../OUTPUT/output.txt", ios::app);
    coutf << "Finalizing TSP..." << endl;
    coutf << "Best travel found: " << endl;
    for (int i = 0; i < _nCities; i++) {
        coutf << _cities(_population(i, 0)).getX() << " " << _cities(_population(i, 0)).getY() << endl;
    }
    coutf << "Loss: " << _loss(0) << endl;
    coutf.close();
}

void saveBestLoss(int gen) {
    ofstream lossFile("../OUTPUT/loss.txt", ios::app);
    if (!lossFile.is_open()) {
        cerr << "Error: Could not open loss.txt file." << endl;
        return;
    }
    
    lossFile << gen << " " << _loss(0) << endl;
    lossFile.close();
}

void saveBestLoss(int gen, int rank) {
    string filename = "../OUTPUT/loss" + to_string(rank) + ".txt";
    ofstream lossFile(filename, ios::app);
    if (!lossFile.is_open()) {
        cerr << "Error: Could not open " << filename << " file." << endl;
        return;
    }
    
    lossFile << gen << " " << _loss(0) << endl;
    lossFile.close();
}


void partialFinalization(int gen) {
    string filename = "../OUTPUT/output" + to_string(gen) + ".txt";
    ofstream coutf(filename);
    coutf << "Partial results after generation " << gen << ": " << endl;
    for (int i = 0; i < _nCities; i++) {
        coutf << _cities(_population(i, 0)).getX() << " " << _cities(_population(i, 0)).getY() << endl;
    }
    coutf << "Loss: " << _loss(0) << endl;
    coutf.close();
}

int getMigrationStep() const {
    return _migrationStep;
}

int getNCities() const {
    return _nCities;
}

arma::mat getCitiesPositions() {
    arma::mat positions(_nCities, _dimension);
    for (int i = 0; i < _nCities; i++) {
        positions.row(i) = _cities(i).getCoordinates().t();
    }
    return positions;
}

void setCitiesPositions(const arma::mat& positions) {
    if (positions.n_rows != _nCities || positions.n_cols != _dimension) {
        cerr << "Error: Dimension mismatch in setCitiesPositions." << endl;
        return;
    }
    for (int i = 0; i < _nCities; i++) {
        arma::vec pos = positions.row(i).t();
        _cities(i).setLocation(pos);
    }
}

arma::Col<int> getbestTravel(int gen) {
    return arma::conv_to<arma::Col<int>>::from(_population.col(0));
}
void set_best_travel(const arma::Col<int>& bestTravel) {
    if (bestTravel.n_elem != _nCities) {
        cerr << "Error: Dimension mismatch in set_best_travel." << endl;
        return;
    }
    for (int i = 0; i < _nCities; i++) {
        _population(i, 0) = bestTravel(i);
    }
}

private:
    Random _rnd;
    const int _dimension = 2; // 2D space
    int _type; //Problem type: 0 for Circle, 1 for scatter
    double _normOrder = 2; // Norm order for distance calculation, default is 2 (Euclidean norm)
    double _power = 2; 
    int _nCities = 10; // Number of cities to be generated
    int _nIndividuals = 100; // Number of individuals in the population
    int _nGenerations = 1000; // Number of generations for the genetic algorithm
    int _migrationStep = 10; // Migration step for the genetic algorithm
    int _rank; // MPI rank
    arma::vec _probMutations; // Probability of mutations
    arma::field<City> _cities; // Field to store cities
    arma::mat _population;
    arma::vec _loss;
};

#endif