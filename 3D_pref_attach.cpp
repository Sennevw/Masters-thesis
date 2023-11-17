#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <array>
#include <vector>
#include <random>
#include <chrono>
#include <map>

using namespace std;
class Particle{
    public:
        int x;
        int y;
        int z;
        int mass;
        Particle(int x, int y, int z, int mass){
            this->x = x;
            this->y = y;
            this->z = z;
            this->mass = mass;
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
    const int size = 50;
    const int iterations = 10;
    unsigned long nparticles = 0;
    const int max_particles = 100000; 
    const int radius = 6;
    const int steps = 40000;

    mt19937 mt(time(0));

    array<array<array<Particle*, size>, size>, size> grid;
    array<Particle*, max_particles> tot_part_list;
    vector<Particle*> used_list;
    vector<float> chance_list;
    vector<Particle*> ptr_list;

    map<int, int> tot_size;
    int count = 0;

    for (int i = 0; i < max_particles; i++){
        Particle* p = new Particle(0, 0 , 0, 1);
        tot_part_list[i] = p;
    }

    for(int iteration = 0; iteration < iterations; iteration++){
        cout << "\0033c \r " << iteration << endl;
        clear_grid<size>(grid);
        for (int i = 0; i < max_particles; i++){
            Particle* p = tot_part_list[i];
            p->x = static_cast<int>(mt() % static_cast<unsigned int>(size));
            p->y = static_cast<int>(mt() % static_cast<unsigned int>(size));
            p->z = static_cast<int>(mt() % static_cast<unsigned int>(size));
            p->mass = 1;
            if (grid[p->x][p->y][p->z] != nullptr){
                grid[p->x][p->y][p->z]->mass ++;
                p->mass = 0;
            }
            else{
                grid[p->x][p->y][p->z] = p;
                used_list.push_back(p);
                nparticles++;
            }
        }
        cout << nparticles << endl;

        for (int i = 0; i < steps; i++){
            int random = mt() % (used_list.size() - 1);
            Particle* p = used_list[random];
            for (int j = -radius; j <= radius; j++){
                for (int k = -radius; k <= radius; k++){
                    for (int l = -radius; l <= radius; l++){
                        if (j == 0 && k == 0 && l == 0) continue;
                        if (grid[(p->x + j) % size][(p->y + k) % size][(p->z + l) % size] == nullptr) continue;
                        if (sqrt(pow(j, 2) + pow(k, 2) + pow(l, 2)) <= radius) {
                            chance_list.push_back((grid[(p->x + j) % size][(p->y + k) % size][(p->z + l) % size]->mass)/(pow(j, 2) + pow(k, 2) + pow(l, 2)));
                            ptr_list.push_back(grid[(p->x + j) % size][(p->y + k) % size][(p->z + l) % size]);
                            
                        }
                        
                    }
                }
            }
            if (chance_list.size() > 0) {
                discrete_distribution<int> dist(chance_list.begin(), chance_list.end());
                Particle* ptr = ptr_list[dist(mt)];
                ptr->mass ++;
                p->mass --;
                if (p->mass == 0) {
                    grid[p->x][p->y][p->z] = nullptr;
                    swap(used_list[random], used_list.back());
                    used_list.pop_back();
                    nparticles --;
                }
                chance_list.clear();
                ptr_list.clear();
            }
        }
        for (int i = 0; i < max_particles; i++){
            ++tot_size[tot_part_list[i]->mass];
        }
        nparticles = 0;
        
    }
    ofstream myfile;
    myfile.open("R" + to_string(radius)+"_PA_"+to_string(max_particles/pow(size, 3))+".txt");
    myfile << iterations << " " << nparticles << " " << radius << endl;
    for (auto const& x : tot_size){
        myfile << x.first << " " << x.second << endl;
    }
    myfile.close();
    return 1;
}