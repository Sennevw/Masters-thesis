#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <array>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
class Particle{
    public:
        int x;
        int y;
        int z;
        int mass;
        Particle* part_of;
        Particle(int x, int y, int z, int mass){
            this->x = x;
            this->y = y;
            this->z = z;
            this->mass = mass;
            this->part_of = nullptr;
        }
};

template <int N>
void clear_grid(array<array<array<Particle*, N>, N>, N>& grid){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                grid[i][j][k] = nullptr;
            }
        }
    }
}

int main(){
    const unsigned long size = 50;
    const int iterations = 500;
    unsigned long nparticles = 0;
    const int max_particles = 100000; 
    const bool same_value = false;
    const int radius = 2;
    const int steps = 50000;

    mt19937 mt(time(0));

    array<array<array<Particle*, size>, size>, size> grid;
    array<Particle*, max_particles> tot_part_list;
    vector<Particle*> used_list;
    vector<float> chance_list;
    vector<Particle*> ptr_list;

    array<int, max_particles> tot_size;

    for (int i = 0; i < max_particles; i++){
        Particle* p = new Particle(0, 0 , 0, 1);
        tot_part_list[i] = p;
    }

    for(int iteration = 0; iteration < iterations; iteration++){
        cout << "\0033c \r " << iteration << endl;
        clear_grid<size>(grid);
        for (int i = 0; i < max_particles; i++){
            Particle* p = tot_part_list[i];
            p->x = mt() % size;
            p->y = mt() % size;
            p->z = mt() % size;
            p->mass = 1;
            if (grid[p->x][p->y][p->z] != nullptr){
                grid[p->x][p->y][p->z]->mass += 1;
                p->mass = 0;
            }
            else{
                grid[p->x][p->y][p->z] = p;
                used_list.push_back(p);
                nparticles++;
            }
        }

        for (int i = 0; i < steps; i++){
            int random = mt() % nparticles;
            Particle* p = used_list[random];
            if (p->part_of != nullptr) {
                p = p->part_of;
            }
            for (int j = -radius; j <= radius; j++){
                for (int k = -radius; k <= radius; k++){
                    for (int l = -radius; l <= radius; l++){
                        if (j !=0 || k != 0 || l != 0){
                            if (grid[(p->x + j) % size][(p->y + k) % size][(p->z + l) % size] != nullptr){
                                if (sqrt(pow(j, 2) + pow(k, 2) + pow(l, 2)) <= radius) {
                                    ptr_list.push_back(grid[(p->x + j) % size][(p->y + k) % size][(p->z + l) % size]);
                                    chance_list.push_back((grid[(p->x + j) % size][(p->y + k) % size][(p->z + l) % size]->mass * p->mass)/(pow(j, 2) + pow(k, 2) + pow(l, 2)));
                                }
                            }
                        }
                    }
                }
            }
            if (ptr_list.size() == 1) {
                Particle* ptr = ptr_list[0];
                ptr->mass += p->mass;
                p->mass = 0;
                grid[p->x][p->y][p->z] = nullptr;
                if (same_value){
                    p->part_of = ptr;
                    used_list[random] = ptr;
                }
                if (!same_value){
                    swap(used_list[random], used_list.back());
                    used_list.pop_back();
                    nparticles--;
                }
                chance_list.clear();
                ptr_list.clear();
            }
            else if (ptr_list.size() != 0) {
                discrete_distribution<int> dist(chance_list.begin(), chance_list.end());
                Particle* ptr = ptr_list[dist(mt)];
                ptr->mass += p->mass;
                p->mass = 0;
                grid[p->x][p->y][p->z] = nullptr;
                if (same_value){
                    p->part_of = ptr;
                    used_list[random] = ptr;
                }
                if (!same_value){
                    swap(used_list[random], used_list.back());
                    used_list.pop_back();
                    nparticles--;
                }
                chance_list.clear();
                ptr_list.clear();
            }
        }
        for (int i = 0; i < max_particles; i++){
            tot_size[tot_part_list[i]->mass] += 1;
            tot_part_list[i]->mass = 0;
        }
    }
    ofstream myfile;
    myfile.open("R" + to_string(radius)+"_PA_"+to_string(max_particles/pow(size, 3))+".txt");
    myfile << iterations << " " << nparticles << " " << radius << endl;
    for (int i = 0; i < max_particles; i++){
        myfile << i << " " << tot_size[i] << endl;
    }
    return 1;
}