#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <array>
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <algorithm>

using namespace std;
template <int N>
void clear_grid(array<array<array<int, N>, N>, N>& grid){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                grid[i][j][k] = 0;
            }
        }
    }
}

void print_vector(vector<float> v){
    for (int i = 0; i < int(v.size()); i++){
        cout << float(v[i]) << " ";
    }
    cout << endl;
}

using namespace std;
int main(){
    const unsigned long size = 50;
    const int nparticles = 100000;
    const float radius = 5;
    int iterations = 500;
    
    // initialize random seed
    mt19937 mt(time(0));

    // generate grid
    array<array<array<int, size>, size>, size> grid;
    // list that will be used store the amount of particles of a certain size
    array<int, nparticles> total_size_list;
    // list that holds the chance of each particle near the last placed particle to be chosen
    vector<float> chance_list;
    // list that holds pointers to all particles near last placed particle
    vector<int*> ptr_list;

    // store random coordinates
    array<int, nparticles> x_list;
    array<int, nparticles> y_list;
    array<int, nparticles> z_list;

    for (int iteration = 0; iteration <iterations; iteration++){
        cout << "\0033c \r " << iteration << endl;
        clear_grid<size>(grid);
        //generate random coordinates
        for (int i = 0; i < nparticles; i++){
            x_list[i] = mt() % size;
            y_list[i] = mt() % size;
            z_list[i] = mt() % size;
        }
        for (int i = 0; i < nparticles; i++){
            // place particle
            grid[x_list[i]][y_list[i]][z_list[i]] += 1;
            // run over box of radius around particle
            for (int j = -radius; j <= radius; j++){
                for (int k = -radius; k <= radius; k++){
                    for (int l = -radius; l <= radius; l++){
                        // do not look at particle itself
                        if (j !=0 || k != 0 || l != 0){
                            // mod is used to make the box periodic
                            // if a particle is found in the box its chance is calculated and added to chance_list, its pointer is added to ptr_list
                            if (grid[(x_list[i]+j) % size][(y_list[i]+k) % size ][(z_list[i]+l) % size] != 0){
                                if (sqrt(pow(j, 2) + pow(k, 2) + pow(l, 2)) <= radius) {
                                    chance_list.push_back(grid[(x_list[i]+j) % size][(y_list[i]+k) % size][(z_list[i]+l) % size] / (pow(j, 2) + pow(l, 2) + pow(k,2)));
                                    ptr_list.push_back(&grid[(x_list[i]+j) % size][(y_list[i]+k) % size][(z_list[i]+l) % size]);
                                }
                            }
                        }
                    }
                }
            }
            // if there are particles in the box, a particle is chosen according to its chance and its size is increased by 1
            if (chance_list.size() == 1){
                int* ptr = ptr_list[0];
                *ptr += grid[x_list[i]][y_list[i]][z_list[i]];
                grid[x_list[i]][y_list[i]][z_list[i]] = 0;
                chance_list.clear();
                ptr_list.clear();
            }
            else if (!chance_list.empty()){
                discrete_distribution<int> dist(chance_list.begin(), chance_list.end());
                int* ptr = ptr_list[dist(mt)];
                *ptr += grid[x_list[i]][y_list[i]][z_list[i]];
                grid[x_list[i]][y_list[i]][z_list[i]] = 0;
                chance_list.clear();
                ptr_list.clear();
            }
        }
        // count the amount of particles of each size
        for (int i = 0; i < int(size); i++){
            for (int j = 0; j < int(size); j++){
                for (int k = 0; k < int(size); k++){
                    total_size_list[grid[i][j][k]] += 1;
                }
            }
        }
    }
    // write data to file
    ofstream myfile;
    myfile.open("r"+to_string(radius) +'_'+to_string(nparticles/pow(size, 3))+".txt");
    myfile << iterations << " " << nparticles << " " << radius << endl;
    for (int i = 0; i < nparticles; i++){
        myfile << i << " " << total_size_list[i] << endl;
    }
    return 1;
}