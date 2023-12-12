#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <array>
#include <map>
#include <iomanip>
#include <list>
#include <sstream>

using namespace std;
template <typename T>
class Point{
    public:
        vector<T> coords;
        int mass;
        static const int Dim = 2;

        Point(){};

        Point(vector<T> coord){
            this->coords = coord;
            this->mass = 1;
        }
        Point(vector<T> coord, int mass){
            this->coords = coord;
            this->mass = mass;
        }
};
template <class T, template<typename> class Point>
float distance(Point<T>& p1, Point<T>& p2, int size){
    int dist = 0;
    for (int i = 0; i < Point<T>::Dim; i++){
        dist += pow( min(fabs(p1.coords[i] - p2.coords[i]), size - fabs(p1.coords[i] - p2.coords[i])), Point<T>::Dim - 1);
    }
    return dist;
}
int main(){
    int iterations = 20;
    int size = 500;
    int nPoints = 25000;
    int nDim = 2;
    int runTime = 25000;
    int measurements = 20;

    bool mode = false;

    vector<Point<int>> points;
    points.reserve(nPoints);

    vector<float> weights;
    weights.reserve(nPoints);
    vector<map<int, int>> mass_list(measurements);
    vector<map<int, int>>::iterator it = mass_list.begin();
    vector<int> indices;


    if (nDim != Point<int>::Dim){throw runtime_error("Dimension of point is not equal to nDim");}

    mt19937 mt(time(0));
    uniform_int_distribution<int> dist(0, size);
    uniform_int_distribution<int> dist2(0, nPoints);


    for(int iteration = 0; iteration < iterations; iteration++) {
        cout << '\r' << string(30, ' ') << '\r' << flush;
        cout << iteration << flush;
        it = mass_list.begin();
        for (int i = 0; i < nPoints; i++){
            vector<int> coord;
            for (int dim = 0; dim < nDim; dim++){
                coord.push_back(static_cast<int>(dist(mt)));
            }
            points.push_back(Point<int>(coord));
        }
        for (int time = 0; time <= runTime; time++){
            int randIdx = static_cast<int>(dist2(mt));
            int loop_count = 0;
            bool to_break = false;
            while(points[randIdx].mass == 0){
                randIdx = static_cast<int>(dist2(mt));
                loop_count++;
                if (loop_count > nPoints) {to_break = true; break;}
            }
            if (to_break) {points.clear(); break;}
            Point<int>& randPoint = points[randIdx];
            for (int i = 0; i < nPoints; i++){
                if (i == randIdx) weights.push_back(0);
                else {weights.push_back( points[i].mass/distance(points[i], randPoint, size));}
            }
            discrete_distribution<int> dist3(weights.begin(), weights.end());
            int idx = dist3(mt);
            switch (mode){
                case true:
                    points[idx].mass += 1;
                    points[randIdx].mass -= 1;
                    break;
                case false:
                    points[idx].mass += randPoint.mass;
                    points[randIdx].mass = 0;
                    break;
            }
            weights.clear();
            if (time != 0  && (time % (runTime/measurements)) == 0) {
                for (int i = 0; i < int(points.size()); i++) {
                    (*it)[points[i].mass] ++;
                }
                it++;
            }
            if (time % 100 == 0) {
                cout << "\r" << string(30, ' ') << flush;
                cout << time/float(runTime) * 100 << "%" << flush;
            }
        }
        points.clear();
    }
    for (int i = 0; i < int(mass_list.size()); i++) {
        string name = "rmax_" + to_string(int(nDim)) + "D_" + to_string(size) + "_" +to_string((i+1) * runTime/measurements) +'_'+to_string(nPoints/pow(size, nDim)) + '_' + to_string(mode) +".txt";
        cout << name << endl;
        ofstream myfile;
        myfile.open("/home/sennevw/Thesis/Pref_attach/" + name);
        myfile << iterations << " " << nPoints << " " << "max" << endl;
        for (const auto& n: mass_list[i])
            myfile << n.first << " " << n.second << endl;
        myfile.close();
        }
}