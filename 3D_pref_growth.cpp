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
#include <map>

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
    const int size = 50;
    const int nparticles = 100000;
    const double radius = 5;
    int iterations = 500;
    
    // initialize random seed
    mt19937 mt(time(0));

    // generate grid
    array<array<array<int, size>, size>, size> grid;
    // list that will be used store the amount of particles of a certain size
    map<int, int> total_size_list;
    // list that holds the chance of each particle near the last placed particle to be chosen
    vector<float> chance_list;
    // list that holds pointers to all particles near last placed particle
    //vector<int*> ptr_list;
    vector<vector<int>> ptr_list2;

    // store random coordinates
    array<int, nparticles> x_list;
    array<int, nparticles> y_list;
    array<int, nparticles> z_list;

    for (int iteration = 0; iteration <iterations; iteration++){
        cout << "\0033c \r " << iteration << endl;
        clear_grid<size>(grid);
        //generate random coordinates
        for (int i = 0; i < nparticles; i++){
            x_list[i] = static_cast<int>(mt() % static_cast<unsigned int>(size));
            y_list[i] = static_cast<int>(mt() % static_cast<unsigned int>(size));
            z_list[i] = static_cast<int>(mt() % static_cast<unsigned int>(size));
        } 
        for (int i = 0; i < nparticles; i++){
            // place particle
            if (grid[x_list[i]][y_list[i]][z_list[i]] != 0) {
                grid[x_list[i]][y_list[i]][z_list[i]] += 1;
                continue;
            }
            // run over box of radius around particle
            for (int j = -radius; j <= radius; j++){
                for (int k = -radius; k <= radius; k++){
                    for (int l = -radius; l <= radius; l++){
                        // do not look at particle itself
                        if (j == 0 && k == 0 && l == 0) continue;
                            // mod is used to make the box periodic
                            // if a particle is found in the box its chance is calculated and added to chance_list, its pointer is added to ptr_list
                        if (grid[(x_list[i]+j) % size][(y_list[i]+k) % size ][(z_list[i]+l) % size] != 0){
                            if (sqrt(pow(j, 2) + pow(k, 2) + pow(l, 2)) <= radius) {
                                chance_list.push_back(grid[(x_list[i]+j) % int(size)][(y_list[i]+k) % int(size)][(z_list[i]+l) % int(size)] / double((pow(j, 2) + pow(l, 2) + pow(k, 2))));
                                ptr_list2.push_back({(x_list[i]+j) % int(size), (y_list[i]+k) % int(size), (z_list[i]+l) % int(size)});
                            }
                        }
                    }
                }
            }
            // if there are particles in the box, a particle is chosen according to its chance and its size is increased by 1
            if (chance_list.empty()) {
                grid[x_list[i]][y_list[i]][z_list[i]] += 1;
            }
            else if (chance_list.size() != 0){
                discrete_distribution<int> dist(chance_list.begin(), chance_list.end());
                vector<int> coord = ptr_list2[dist(mt)];
                grid[coord[0]][coord[1]][coord[2]] += 1;
                chance_list.clear();
                ptr_list2.clear();
            }
        }
        // count the amount of particles of each size
        for (int i = 0; i < int(size); i++){
            for (int j = 0; j < int(size); j++){
                for (int k = 0; k < int(size); k++){
                    ++total_size_list[grid[i][j][k]];
                }
            }
        }
    }   
    // write data to file
    cout << "r"+to_string(radius) +'_'+to_string(nparticles/pow(size, 3))+".txt" << endl;

    ofstream myfile;
    myfile.open("r"+to_string(radius) +'_'+to_string(nparticles/pow(size, 3))+".txt");
    myfile << iterations << " " << nparticles << " " << radius << endl;
    for (auto const& x : total_size_list){
        myfile << x.first << " " << x.second << endl;
    }
    myfile.close();
    return 0;
}