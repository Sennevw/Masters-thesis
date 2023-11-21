#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <map>

using namespace std;
struct Particle {
    int id;
    double inv_mass;
    Particle(int id, double inv_mass) : id(id), inv_mass(inv_mass) {}
};

int main() {
    int numSteps = 1000;
    int iterations = 5000;

    vector<Particle> particles;
    particles.reserve(numSteps + 1);

    map<int, int> particleMap;

    bool mode = false;
    double lambda = 1;
    uniform_real_distribution<double> dist(0.0, 1.0);
    double split;
    
    random_device rd;
    mt19937 gen(rd());
    
    for (int iteration = 0; iteration < iterations; iteration++){
        particles.push_back(Particle(1, 1)); 
        Particle& selectedParticle = particles[0];
        for (int step = 0; step < numSteps; step++){
            int index = static_cast<int>(gen() % static_cast<unsigned int>((particles.size()))); 
            selectedParticle = particles[index];
            particles.push_back(Particle(int(particles.size()) + 1, selectedParticle.inv_mass * 0.5));
            selectedParticle.inv_mass *= 0.5;
        }
        cout << "\r" << string(30, ' ') << flush;
        cout << iteration/float(iterations) * 100 << "%" << flush;
        for (auto& particle : particles) {
            particleMap[particle.inv_mass] += 1;
        }
        particles.clear();
        }
    
    ofstream myfile;
    myfile.open("alternative" + to_string(numSteps) + '_' + to_string(mode) + ".txt");
    for (auto& [key, value] : particleMap) {
        myfile << key << " " << value << endl;
    }
    myfile.close();
    return 0;
}
