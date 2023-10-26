#include <cmath>
#include "boost/multi_array.hpp"
#include <iostream>
#include <array>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <random>
#include <algorithm>



using namespace std;
template <int N>
void clear_field(array<array<array<int, N>, N>, N>& field, int background){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                field[i][j][k] = background;
            }
        }
    }
}


int main(){
    const int size = 100;
    const int nparticles = 10000;
    int iterations = 1;
    int background = 10;

    array<array<array<int, size>, size>, size> densfield;
    array<array<array<int, size>, size>, size> densfieldnew;
    array<int, nparticles> total_size_list;

    //fill densfield with constant background
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            for (int k = 0; k < size; k++){
                densfield[i][j][k] = background;
            }
        }
    }

    array<int, nparticles> x_list;
    array<int, nparticles> y_list;
    array<int, nparticles> z_list;

    int denshere;
    double initvelo;
    int cellsperdim;
    double densmax;

    array<double, 3> velo;
    array<double, 3> dist_rad;
    uniform_real_distribution<double> dist(0.0, 1.0);

    for (int iteration = 0; iteration < iterations; iteration++){
        
        mt19937 mt(random_device{}());
        for (int i = 0; i < nparticles; i++){
            x_list[i] = mt() % size;
            y_list[i] = mt() % size;
            z_list[i] = mt() % size;
        }

        for (int i = 0; i < nparticles; i++){
            densfield[x_list[i]][y_list[i]][z_list[i]] += 1;
            for (int x = 0; x < size; x++){
                for (int y = 0; y < size; y++){
                    for (int z = 0; z < size; z++){
                        denshere = densfield[x][y][z];
                        for (int constituent = 0; constituent < denshere; constituent++){
                            for (int j = 0; j < 3; j++){
                                velo[j] = initvelo*cellsperdim*(1-min(denshere/densmax, 1.0));
                            }
                        }   
                    }
                }
            }
        }
    }
}