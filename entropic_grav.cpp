#include <cmath>
#include "boost/multi_array.hpp"
#include <iostream>
#include <array>



using namespace std;

// wat zijn particlefield en particleparts

void main() {
    const size_t ndim = 200;
    const int len = 100;
    int iterations = 1;
    double particleparts = 0.0; // ??????????

    array<int, ndim> dim;
    for (size_t i = 0; i < ndim; i++) {
        dim[i] = len;
    }
    boost::multi_array<double, ndim> densfield(dim);
    boost::multi_array<double, ndim> densfieldnew(dim);
    boost::multi_array<double, ndim> particlefield(dim);
    for (int i = 0; i < iterations; i++) {
        for (int j = 0; j < pow(len,ndim); j++){
            densfield(dim) = 0;
            // run over all elements in densfield
            for (int k = ndim - 1; k > -1; k--){   
            }
        }
        while(dim[0] <len){
            int denshere = 0;
            int densout = 0;
            if (particlefield(dim) == 0){
                densout = denshere - particleparts;

            }

        }
    }
}






