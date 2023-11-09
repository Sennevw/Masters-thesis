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
void clear_grid(array<array<int, N>, N>& grid){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            grid[i][j] = 0;
        }
    }
}


using namespace std;
int main(){
    const unsigned long size = 50;
    const int nparticles = 2000;
    const float radius = 5;
    int iterations = 500;
    
    // initialize random seed
    mt19937 mt(time(0));

    // generate grid
    array<array<int, size>, size> grid;
    // list that will be used store the amount of particles of a certain size
    
    // list that holds the chance of each particle near the last placed particle to be chosen
    vector<float> chance_list;
    // list that holds pointers to all particles near last placed particle
    //vector<int*> ptr_list;
    vector<vector<int>> ptr_list2;
    // list that holds the coordinates of all particles
    vector<int> x_list(nparticles);
    vector<int> y_list(nparticles);
    array<int, nparticles> total_size_list;
    fill(total_size_list.begin(), total_size_list.end(), 0);


    // store random coordinates
    for (int iteration = 0; iteration < iterations; iteration++){
        cout << iteration << endl;
        clear_grid<size>(grid);
        //generate random coordinates
        for (int i = 0; i < nparticles; i++){
            x_list[i] = mt() % size;
            y_list[i] = mt() % size;
        } 

        for (int i = 0; i < nparticles; i++){
            // place particle
            if (grid[x_list[i]][y_list[i]] != 0) {
                grid[x_list[i]][y_list[i]] += 1;
                continue;
            }
            // run over box of radius around particle
            for (int j = -radius; j <= radius; j++){
                for (int k = -radius; k <= radius; k++){
                    // do not look at particle itself
                        // mod is used to make the box periodic
                        // if a particle is found in the box its chance is calculated and added to chance_list, its pointer is added to ptr_list
                    if (grid[(x_list[i]+j) % size][(y_list[i]+k) % size ] != 0){
                        if (sqrt(pow(j, 2) + pow(k, 2)) <= radius) {
                            chance_list.push_back(grid[(x_list[i]+j) % int(size)][(y_list[i]+k) % int(size)] / (pow(j, 2) + pow(k, 2) ));
                            ptr_list2.push_back({(x_list[i]+j) % int(size), (y_list[i]+k) % int(size)});
                        }
                    }
                }
            }
            // if there are particles in the box, a particle is chosen according to its chance and its size is increased by 1
            if (chance_list.empty()) {
                grid[x_list[i]][y_list[i]] += 1;
            }
            else if (chance_list.size() != 0){
                discrete_distribution<int> dist(chance_list.begin(), chance_list.end());
                vector<int> coord = ptr_list2[dist(mt)];
                grid[coord[0]][coord[1]] += 1;
                chance_list.clear();
                ptr_list2.clear();
            }
        }
        // count the amount of particles of each size
        for (int i = 0; i < int(size); i++){
            for (int j = 0; j < int(size); j++){
                total_size_list[grid[i][j]] += 1;
            }
        }

    }
    // write data to file
    
    std::cout << "2D_r"+to_string(radius) +'_'+to_string(nparticles/pow(size, 2))+".txt" << endl;
    ofstream myfile;
    myfile.open("2D_r"+to_string(radius) +'_'+to_string(nparticles/pow(size, 2))+".txt");
    myfile << iterations << " " << nparticles << " " << radius << endl;
    for (int i = 0; i < nparticles; i++){
        myfile << i << " " << total_size_list[i] << endl;
    }
    myfile.close();
    return 0;
};