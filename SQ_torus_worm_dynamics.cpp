//The code is already written in python. Translate that to c++ 
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <array> 
#include <complex>
#include <cmath>
#include <algorithm> // for std::sort
//#include <Eigen> 
#include <random> 
#include <chrono>
// #include <python.h>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <future>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;
using namespace std::chrono;
//using namespace Eigen;

// //First define the square lattice with random couplings

//call the data stored in csv files
vector<vector<int>> readMatrixFromCSV(const string& filename) {
    vector<vector<int>> matrix;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Unable to open file" << endl;
        return matrix;
    }

    string line;
    while (getline(file, line)) {
        vector<int> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, ',')) {
            row.push_back(stod(value));
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

vector<vector<double>> readMatrixFromCSV_doub(const string& filename) {
    vector<vector<double>> matrix;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Unable to open file" << endl;
        return matrix;
    }

    string line;
    while (getline(file, line)) {
        vector<double> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, ',')) {
            row.push_back(stod(value));
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}



vector<int> ad_x(int x, int N) {
    if (x == 1) {
        return {N, 2};
    } else if (x == N) {
        return {N - 1, 1};
    } else {
        return {x - 1, x + 1};
    }
}

vector<int> ad_y(int y, int M) {
    if (y == 1) {
        return {M, 2};
    } else if (y == M) {
        return {M - 1, 1};
    } else {
        return {y - 1, y + 1};
    }
}
int Enum(int x, int y, int k, int l, int N, int M) {
    int v;
    v = 6*(N*(y-1)+(x-1))+3*(k-1)+l;
    return v;
}

vector<int> Rev_Enum(int c, int N, int M) {
    int x, y, k, l;
    int m1 = std::ceil(static_cast<double>(c ) / (6 * N));
    int n1 = std::ceil(static_cast<double>(c - (6 * (m1 - 1)) * N) / 6);
    x = n1;
    y = m1;
    k = std::ceil(static_cast<double>(c - (6 * N * (m1 - 1) + 6 * (n1 - 1))) / 3);
    l = c - (6 * N * (m1 - 1) + 6 * (n1 - 1)) - 3 * (k - 1);
    return {x, y, k, l};
}

float Nishimori(float b){
    float p;
    float a=1+exp(-2*b);
    p=pow(a,-1);
    
  return p;
}

//This saves the Neigbours matrix by calling it from the .csv file
std::vector<std::vector<int>> Neighbours;
void initializeNeighbours() {
    std::string filename = "NEIG.csv";
    std::vector<std::vector<int>> matrix = readMatrixFromCSV(filename);

    int v = matrix.size();
    int v1 = matrix[0].size();
    Neighbours.resize(v, std::vector<int>(v1));

    // Convert the matrix elements to integers and store them in Neighbours
    for (int i = 0; i < v; ++i) {
        for (int j = 0; j < v1; ++j) {
            Neighbours[i][j] = static_cast<int>(matrix[i][j]);
        }
    }
}

// The Neigh function finds the neigbours of some point
std::vector<int> Neigh(int c) {
    // Ensure the Neighbours matrix is initialized
    if (Neighbours.empty()) {
        std::cerr << "Neighbours matrix has not been initialized!" << std::endl;
        return {};
    }

    // Check for a valid index
    if (c - 1 < 0 || c - 1 >= Neighbours.size()) {
        std::cerr << "Invalid index for Neighbours: " << c - 1 << std::endl;
        return {};
    }

    // Return the row corresponding to the index
    return Neighbours[c - 1];
}

//Initial matching
std::vector<std::vector<int>> Matching0;
void initializeMatchings() {
    std::string filename = "Mat.csv";
    std::vector<std::vector<int>> matrix = readMatrixFromCSV(filename);

    int v = matrix.size(); //Can be simply written, this will make it faster
    int v1 = matrix[0].size();
    Matching0.resize(v, std::vector<int>(v1));

    // Convert the matrix elements to integers and store them in Neighbours
    for (int i = 0; i < v; ++i) {
        for (int j = 0; j < v1; ++j) {
            Matching0[i][j] = static_cast<int>(matrix[i][j]);
        }
    }
}

//Detailed balance matrix
std::vector<std::vector<double>> DBM;
void initializeDBM() {
    std::string filename = "Det_balance.csv";
    std::vector<std::vector<double>> matrix = readMatrixFromCSV_doub(filename);

    int v = matrix.size();
    int v1 = matrix[0].size();
    DBM.resize(v, std::vector<double>(v1));

    // Convert the matrix elements to integers and store them in Neighbours
    for (int i = 0; i < v; ++i) {
        for (int j = 0; j < v1; ++j) {
            DBM[i][j] = static_cast<double>(matrix[i][j]);
        }
    }
}

double randomDouble(){
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

// Dimer move function
//x2 is the pivot site
int Dimer_move(int x1, int x2) {

    // Adjust for 1-based to 0-based index (C++ uses 0-based indexing)
    std::vector<double> A1 = DBM[x2 - 1];        // Extract the row of matrix A for the pivot x2
    std::vector<int> NB = Neighbours[x2 - 1];  // Neighbours of the pivot x2

    // Adjust the neighbour list, assuming that NB has at least 4 elements as in Python code
    std::vector<int> NB_ = {NB[1], NB[2], NB[3]}; // Keep only three neighbours

    // Find the index of x1 in NB_
    auto it = std::find(NB_.begin(), NB_.end(), x1);
    if (it == NB_.end()) {
        std::cerr << "x1 not found in the neighbours of x2!" << std::endl;
        return -1;  // Return an error value
    }
    int ind = std::distance(NB_.begin(), it);  // Index of x1 in NB_

    // Calculate transition probabilities
    std::vector<double> P = {
        A1[ind * 3] / A1[9 + ind],
        A1[ind * 3 + 1] / A1[9 + ind],
        A1[ind * 3 + 2] / A1[9 + ind]
    };

    // Generate a random number between 0 and 1
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Discard the first number (optional warm-up step)
    distribution(generator);

    // Generate the next random number
    double r = distribution(generator);

    // Perform the transition based on the probabilities P
    if (r <= P[0]) {
        return NB_[0];  // Return the first neighbour
    } else if (r <= P[0] + P[1]) {
        return NB_[1];  // Return the second neighbour
    } else {
        return NB_[2];  // Return the third neighbour
    }
}


int Dimer_identifier(int a, const vector<vector<int>>& Matching_) {
    for (const auto& s : Matching_) {
        if (s[0] == a) {
            return s[1];
            break;
        } else if (s[1] == a) {
            return s[0];
            break;
        }
    }
    return 0; // Return -1 if no match is found
}

void printList(const std::vector<int>& List) {
    for (const auto& d : List) {
        std::cout << d << std::endl;
    }
}


int main() {

    std::vector<int> Distance;

    ofstream outputFile("Len_dist_bc_10000_L101_Square_p100.csv", ios_base::app);
    if (!outputFile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return 1;
    }
    
    initializeNeighbours();
    initializeDBM();
    initializeMatchings();

    double Ensemble = 500;
    double runs = 10000;
    int N = 101, M = 101;

    auto start = chrono::high_resolution_clock::now();
    int x10 = 2000, x20 = Dimer_identifier(2000, Matching0); // Updated initial Matching_0, starting location does not matter  
    vector<int> rev10 = Rev_Enum(x10, N, M);
    int X10 = rev10[0], Y10 = rev10[1];

    double Cn=0;

    for (int en = 1; en <= Ensemble; en++) {
    	vector<vector<int>> Matching = Matching0;  // Create a copy of Matching0
        int length = 0;
        int x1 = x10, x2 = x20;
        int cx=0, cy=0;
        
    for (int i = 1; i <= runs; i++) {

        vector<int> rev1_p = Rev_Enum(x1, N, M);
        int X1_p = rev1_p[0], Y1_p = rev1_p[1];

        int c1 = Dimer_move(x1, x2); // New dimer worm head, first bounce may give the same thing back

        if (c1 == x1) {
            continue; // x2 remains unchanged
        }
        else{

            length += 1;

            Matching.erase(std::remove_if(Matching.begin(), Matching.end(), //Remove the old dimer 
                [&](const std::vector<int>& pair) {
                    return (pair == std::vector<int>{x1, x2} || pair == std::vector<int>{x2, x1});
                }), Matching.end());
                
            Matching.push_back({x2, c1}); //Create the new dimer 
            Matching.push_back({c1, x2});

            x1=c1;

            vector<int> rev1 = Rev_Enum(x1, N, M);
            int X1 = rev1[0], Y1 = rev1[1];

            if (X1 == ad_x(X1_p, N)[1]) { // Counting winding around x
                cx += 1;
            } else if (X1 == ad_x(X1_p, N)[0]) {
                cx -= 1;
            }

            if (Y1 == ad_y(Y1_p, M)[1]) { // Counting winding around y
                cy += 1;
            } else if (Y1 == ad_y(Y1_p, M)[0]) {
                cy -= 1;
            }

            if(x1!=x10){
                int c2 = Dimer_identifier(c1, Matching); // New dimer pivot
                x2=c2;
            }else{
                break;
                
            }

        }

    }
        //Fix the issue of odd length and +-1 counting for cx,cy in some cases. These arises only in (0,0) Homology class.
        double Wx=cx/N, Wy=cy/M;

        if(x1==x10){  
            if ((static_cast<int>(Wx) % 2 == 0) && (static_cast<int>(Wy) % 2 == 0)) {
                //cout<< length<< "," << x1 << "," << x2 << ","<< Wx << "," <<Wy<< endl;
                outputFile << length << "\t"; //Loops in the even, even Homology class
            }else{
                outputFile << 0 << "\t"; //Loops in the wrong Homology classes
            }
        }else{
            outputFile << -1 << "\t"; //The wormhead does not come to the intial defect points
        }

    }

    outputFile << endl;
    outputFile.close();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    //cout<< "Time taken for Sturm sequence method: " << duration.count() << " seconds" << endl;

    return 0;
}
