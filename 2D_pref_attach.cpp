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

int main() {
    const unsigned long size = 50;
    const int nparticles = 100000;
    const int radius = 4;
    int iterations = 1000;

    mt19937 mt(time(0));

    array<array<int, size>, size> grid;
    array<int, nparticles> total_size_list;
    vector<float> chance_list;
    vector<int*> ptr_list;

    array<int, nparticles> x_list;
    array<int, nparticles> y_list;
    
    for (int iteration = 0; iteration < iterations; iteration++){
        cout << iteration << endl;
        clear_grid<size>(grid);
        for (int i = 0; i < nparticles; i++){
            x_list[i] = mt() % size;
            y_list[i] = mt() % size;
        }
        
    }
}