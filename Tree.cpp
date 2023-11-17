#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <array>
#include <iostream>
#include <map>
#include <iomanip>
#include <list>
using namespace std;

struct Node {
    Node* left;
    Node* right;
    int idx;
    int axis;
    Node(): idx(-1), axis(-1){this->left = this->right =nullptr;};
};

template <class T, template<typename> class Point>
struct Neighbour{
    Point<T>& point;
    double distance;
    Neighbour(Point<T>& p, double d): point(p), distance(d){};
};

template<class T, template<typename> class Point>
class Tree{
    public:
        Node* root_;
        vector<Point<T>>* points_;
        float limit;

        Tree(int lim): root_(nullptr), points_(nullptr), limit(lim){};

        Tree(vector<Point<T>>& points, int lim): root_(nullptr), limit(lim){
            if (!points.empty()) {
                points_ = &points;
                build(points);
            }
            else {
                points_ = &points;
            }
        }

        ~Tree() { clear(); }
        
        void build(vector<Point<T>>& points) {
            clear();
            points_ = &points;
            vector<int> indices(points.size());
            iota(indices.begin(), indices.end(), 0);
            root_ = build_recurs(indices.data(), (int)points.size(), 0);
        }

        void clear() {
            clear_recurs(root_);
            root_ = nullptr;
            points_ = nullptr;
        }

        void clear_recurs(Node* node){
            if (node == nullptr) return;
            clear_recurs(node->left);
            clear_recurs(node->right);
            delete node;
        }

        Node* build_recurs(int* indices, int npoints, int depth){
            if (npoints <=0) return nullptr;
            const int axis = depth % Point<T>::Dim;
            const int middle = npoints / 2;

            nth_element(indices, indices + middle, indices + npoints, [&](int lhs, int rhs) {
                return (*points_)[lhs].coords[axis] < (*points_)[rhs].coords[axis];
            });

            Node* node = new Node();
            node->idx = indices[middle];
            node->axis = axis;

            node->left = build_recurs(indices, middle, depth + 1);
            node->right = build_recurs(indices + middle + 1, npoints - middle - 1, depth + 1);

            return node;
        }

        static double distance(const Point<T>& a, const Point<T>& b){
            double dist = 0;
            for (int i = 0; i < Point<T>::Dim; ++i){
                dist += (a.coords[i] - b.coords[i]) * (a.coords[i] - b.coords[i]);
            }
            return sqrt(dist);
        }

        static double pow_distance(const Point<T>& a, const Point<T>& b){
            double dist = 0;
            for (int i = 0; i < Point<T>::Dim; ++i){
                dist += (a.coords[i] - b.coords[i]) * (a.coords[i] - b.coords[i]);
            }
            return dist;
        }

        double circ_distance(const Point<T>& a, const Point<T>& b) {
            double dist = 0;
            for (int i = 0; i < Point<T>::Dim; i++){
                dist += pow( min(fabs(a.coords[i] - b.coords[i]), limit - fabs(a.coords[i] - b.coords[i])), 2);
            }
            return sqrt(dist);
        }

        bool same_coord(const Point<T>& a, const Point<T>& b){
            for (int i = 0; i < Point<T>::Dim; i++){
                if (a.coords[i] != b.coords[i]) return false;
            }
            return true;
        }

        bool are_same(Node* trial, Point<T>& center){
            if (same_coord((*points_)[trial->idx], center)) {

                return true;
            }
            return false;
        }

        bool insert(Point<T>& point){
            if (root_ == nullptr) {
                root_ = new Node();
                root_->idx = 0;
                root_->axis = 0;
                return true;
            }
            else {
                return insert_recurs(root_, point, 0);
            }
        }

        bool insert_recurs(Node* node, Point<T>& insert, int depth){
            const int axis = depth % Point<T>::Dim;
            if (are_same(node, insert)){
                (*points_)[node->idx].mass += 1;
                insert.mass = 0;
                (*points_).pop_back();
                return false;
            }
                
            if (insert.coords[axis] < (*points_)[node->idx].coords[axis]){
                if (node->left == nullptr){
                    node->left = new Node();
                    node->left->idx = points_->size() - 1;
                    node->left->axis = (depth + 1) % Point<T>::Dim;
                    return true;
                }
                else {
                    return insert_recurs(node->left, insert, depth + 1);
                }
            }
            else {
                if (node->right == nullptr){
                    node->right = new Node();
                    node->right->idx = points_->size() - 1;
                    node->right->axis = (depth + 1) % Point<T>::Dim;
                    return true;
                }
                else {
                    return insert_recurs(node->right, insert ,depth + 1);
                }
            }
            return true;
        }

        bool radius_search(Point<T>& center, double radius, vector<int>& indices, bool duplicates){
            bool flip = true;
            radius_search_recurs(root_, center, radius, indices, flip, duplicates);
            Point<T> virtual_center = center;
            vector<vector<T>> coord_list;
            if (!flip) return flip;
            for (int i = 0; i < Point<T>::Dim; i++){
                if (limit - center.coords[i] <= radius){
                    coord_list.push_back(center.coords);
                    coord_list.back()[i] = center.coords[i] - limit;
                    }
                else if (center.coords[i] <= radius) {
                    coord_list.push_back(center.coords);
                    coord_list.back()[i] = center.coords[i] + limit;
                }
            }
            if (coord_list.size() == 1) {
                virtual_center.coords = coord_list[0];
                radius_search_recurs(root_, virtual_center, radius, indices, flip, false);
            }
            else if (coord_list.size() > 1) {
                for (int i = 0; i < Point<T>::Dim; i++) {
                    for (int j = 0; j < int(coord_list.size()); j++){
                        for (int k = 0; k < Point<T>::Dim; k++){
                            if (coord_list[j][k] > 0 && coord_list[j][k] < limit && limit - coord_list[j][k] <= radius){
                                coord_list.push_back(coord_list[j]);
                                coord_list.back()[k] = coord_list[j][k] - limit;
                            }
                            else if (coord_list[j][k] > 0 && coord_list[j][k] <= radius) {
                                coord_list.push_back(coord_list[j]);
                                coord_list.back()[k] = coord_list[j][k] + limit;
                            }
                        }
                    }
                }
                sort(coord_list.begin(), coord_list.end());
                coord_list.erase(unique(coord_list.begin(), coord_list.end()), coord_list.end());
                for (int i = 0; i < int(coord_list.size()); i++){
                    virtual_center.coords = coord_list[i];
                    radius_search_recurs(root_, virtual_center, radius, indices, flip, false);
                }
            }
            return flip;
        };
           

        void circ_radius_search_recurs(Node* trial, Point<T>& center, double radius, vector<int>& indices, bool& flip, bool duplicates){
            if (trial == nullptr) return;
            if (!flip) return;
            if (duplicates && are_same(trial, center)) {
                flip = false;
                (*points_)[trial->idx].mass += 1;
                center.mass = 0;
                (*points_).pop_back();
                }
            if (!flip) return;
            const Point<T>& trial_p = (*points_)[trial->idx];
            if (circ_distance(trial_p, center) <= radius && !are_same(trial_p, center)){
                indices.push_back(trial->idx);
            }
            const int dir = center.coords[trial->axis] < trial_p.coords[trial->axis] ? 0 : 1;
            if (!dir) {
                radius_search_recurs(trial->left, center, radius, indices, flip, duplicates);
            }
            else {
                radius_search_recurs(trial->right, center, radius, indices, flip, duplicates);
            }
            
            const int axis = trial->axis;
            const double diff = fabs(center.coords[axis] - trial_p.coords[axis]);
            if ((diff < radius) || ((limit - diff) < radius)) {
                radius_search_recurs(dir ? trial->left : trial->right, center, radius, indices, flip, duplicates);
            }
            return;
        }

        void radius_search_recurs(Node* trial, Point<T>& center, double radius, vector<int>& indices, bool& flip, bool duplicates){
            if (trial == nullptr) return;
            if (!flip) return;
            if (duplicates && are_same(trial, center)){ 
                flip = false;
                (*points_)[trial->idx].mass += 1;
                center.mass = 0;
                (*points_).pop_back();
                }
            if (!flip) return;
            const Point<T>& trial_p = (*points_)[trial->idx];
            if (trial_p.mass != 0 && distance(trial_p, center) <= radius && !are_same(trial, center)){
                indices.push_back(trial->idx);
            }
            const int dir = center.coords[trial->axis] < trial_p.coords[trial->axis] ? 0 : 1;
            radius_search_recurs(dir ? trial->right : trial->left, center, radius, indices, flip, duplicates);
            const int axis = trial->axis;
            const double diff = fabs(center.coords[axis] - trial_p.coords[axis]);
            if (diff < radius) {
                radius_search_recurs(dir ? trial->left : trial->right, center, radius, indices, flip, duplicates);
            }
            return; 
        }

        void iterate_tree(Node* node, map<int, int>& mass_list){
            if (node == nullptr) return;
            ++mass_list[(*points_)[node->idx].mass];
            iterate_tree(node->left, mass_list);
            iterate_tree(node->right, mass_list);
        }

        void find_neighbours(vector<int>& indices, vector<vector<Neighbour<T, Point>>>& neighbour_list, float radius){
            for (int i = 0; i < (*points_).size(); i++){
                radius_search((*points_)[i], radius, indices, false);
                for (int j = 0; j < indices.size(); j++){
                    neighbour_list[i].push_back(Neighbour<T, Point>((*points_)[indices[j]], circ_distance((*points_)[i], (*points_)[indices[j]])));
                }
                cout << "\r" << string(30, ' ') << flush;
                cout << i/float((*points_).size()) * 100 << "%" << flush;
                indices.clear();
            }
            cout << "\r" << string(30, ' ') << flush;
            cout << "\r" <<"Neighbours found" << endl;
        }

        void gen_rand_point_growth(vector<int>& indices, vector<float>& weights, Point<T>& point, mt19937& mt) {
            for (int j = 0; j < int(indices.size()); j++) {
                    weights.push_back((*points_)[indices[j]].mass/ pow(circ_distance((*points_)[indices[j]], point), 2));
                } 
                discrete_distribution<int> dist(weights.begin(), weights.end());
                (*points_)[indices[dist(mt)]].mass ++;
                (*points_).pop_back();
                
                weights.clear();
                indices.clear();
        }

        void gen_rand_point_attach(vector<int>& indices, vector<float>& weights, Point<T>& point, mt19937& mt, bool mode) {
            for (int j = 0; j < int(indices.size()); j++) {
                    weights.push_back((*points_)[indices[j]].mass/ pow(circ_distance((*points_)[indices[j]], point), 2));
                } 
                discrete_distribution<int> dist(weights.begin(), weights.end());
                
                if (mode) {
                    point.mass--;
                    (*points_)[indices[dist(mt)]].mass ++;}
                else {
                    (*points_)[indices[dist(mt)]].mass += point.mass;
                    point.mass = 0;
                }
                
                weights.clear();
                indices.clear();
        }

        void output(string name, map<int, int>& mass_list, int iterations, int npoints, float radius, int size, int dim, bool mode) {
            cout << name << endl;
            ofstream myfile;
            if (mode) myfile.open("/home/sennevw/Thesis/Pref_attach/" + name);
            else myfile.open("/home/sennevw/Thesis/Pref_growth/" + name);
            myfile << iterations << " " << npoints << " " << radius << endl;
            for (const auto& n: mass_list)
                myfile << n.first << " " << n.second << endl;
            myfile.close();
        }
};



template <typename T>
class Point{
    public:
        vector<T> coords;
        int mass;
        static const int Dim = 3;

        Point(){};

        Point(vector<T> coord){
            this->coords = coord;
            this->mass = 1;
        }        
};

template <class T, template<typename> class Point>
void create_random_float(vector<vector<T>>& coord_list, int npoints, mt19937& mt, int size){
    uniform_real_distribution<float> dist_real(0, size);
    for (int i = 0; i < npoints; i++) {
        for (int j = 0; j < Point<T>::Dim; j++){
            coord_list[i][j] = dist_real(mt);   
        }
    }
}

template <class T, template<typename> class Point>
void create_random_int(vector<vector<T>>& coord_list, int npoints, mt19937& mt, int size){
            for (int i = 0; i < npoints; i++) {
                for (int j = 0; j < Point<T>::Dim; j++){
                    coord_list[i][j] = static_cast<int>(mt() % static_cast<unsigned int>(size));
                }
            }
        }

template <class T, template<typename> class Point>
void gen_init_cond(vector<Point<T>>& points, vector<vector<int>>& coord_list, int npoints, mt19937& mt, int size) {
            create_random_int<T, Point>(coord_list, npoints, mt, size);
            sort(coord_list.begin(), coord_list.end());
            vector<vector<int>>::iterator it1 = coord_list.begin();
            vector<vector<int>>::iterator it2 = coord_list.begin() + 1;
            for (int i = 0; i < npoints; i++){
                while(it2 != coord_list.end() - 1 && *it2 == *it1) it2++;
                points.push_back(Point<T>(*it1));
                points.back().mass += distance(it1, it2)-1;
                
                if (it2 != coord_list.end() - 1) {
                    it1 = it2;
                    it2 = it2 + 1;
                }
                else if (it2 == coord_list.end() - 1) {
                    break;
                }
            }
            cout << "Points found" << endl;
        }

template <typename T>
void pref_growth(int npoints, int size, int dim, float radius, int iterations, mt19937 &mt){
    vector<Point<T>> points;
    points.reserve(npoints);
    vector<int> indices;
    vector<float> weights;
    Tree<T, Point> tree(points, float(size));
    map<int, int> mass_list;
    vector<vector<T>> coord_list(npoints, vector<T>(dim));

    for (int iteration = 0; iteration < iterations; iteration++) {
        create_random_int<T, Point>(coord_list, npoints, mt, size);  
        cout << '\r' << string(30, ' ') << '\r' << flush;     
        cout << iteration << flush;
        for (int i = 0; i < npoints; i++) {
            points.push_back(Point<T>(coord_list[i]));
            if (tree.radius_search(points.back(), radius, indices, true)) {
                if (indices.size() > 0) tree.gen_rand_point_growth(indices, weights, points.back(), mt);
                else tree.insert(points.back());
            }
        }
        tree.iterate_tree(tree.root_, mass_list);
        tree.clear();
        points.clear();
        tree.points_ = &points;
    }
    string name = to_string(dim) + "D_tree_r"+to_string(radius) +'_'+to_string(npoints/pow(size, dim))+".txt";
    tree.output(name, mass_list, iterations, npoints, radius, size, dim, false);
    
};

template <typename T>
void pref_attach(int npoints, int size, int dim, float radius, int iterations, mt19937 &mt, int time, bool mode){
    
    int measurements = 10;
    vector<vector<int>> coord_list(npoints, vector<int>(Point<T>::Dim));
    vector<int> indices;
    vector<float> weights;
    vector<map<int, int>> mass_list(measurements);
    vector<map<int, int>>::iterator it = mass_list.begin();

    vector<Point<T>> points;
    points.reserve(npoints);

    Tree<int, Point> tree(points, float(size));

    int counter = 0;
    int rand_idx = static_cast<int>(mt() % static_cast<unsigned int>(points.size() - 1));
    
    for (int iteration = 0; iteration < iterations; iteration++) {
        it = mass_list.begin();
        cout << '\r' << string(30, ' ') << '\r' << flush;
        cout << iteration << flush;
        create_random_int<T, Point>(coord_list, npoints, mt, size);
        for (int i = 0; i < npoints; i++){
            points.push_back(Point<T>(coord_list[i]));
            tree.insert(points.back());
        }
        cout << points.size() << endl;
        for (int steps = 0; steps <= time; steps++) {
            rand_idx = static_cast<int>(mt() % static_cast<unsigned int>(points.size() - 1));
            while(points[rand_idx].mass == 0) rand_idx = static_cast<int>(mt() % static_cast<unsigned int>(points.size() - 1));
            Point<T>& point = points[rand_idx];
            tree.radius_search(point, radius, indices, false);
            if (indices.size() > 0) tree.gen_rand_point_attach(indices, weights, point, mt, mode);
            if (counter == 1000) {
                tree.clear();
                for (int i = 0; i < int(points.size()); i++) {
                    if (points[i].mass == 0) {
                        swap(points[i], points.back());
                        points.pop_back();
                    }
                }
                tree.build(points);
                counter = 0;
            }
            if (steps != 0  && (steps % (time/10)) == 0) {
                tree.iterate_tree(tree.root_, *it);
                it++;
            }
            cout << "\r" << string(30, ' ') << flush;
            cout << "\r" << steps/float(time) * 100 << "%" << flush;
        }
        indices.clear();
        weights.clear();
        points.clear();
        tree.clear();
        tree.points_ = &points;
        
    }
    for (int i = 0; i < int(mass_list.size()); i++) {
        string name = "r" + to_string(int(radius)) + "_" + to_string(int(dim)) + "D_tree_" +to_string(i * time/10) +'_'+to_string(npoints/pow(size, dim)) + '_' + to_string(mode) +".txt";
        tree.output(name, mass_list[i], iterations, npoints, radius, size, dim, true);
    }
    
};

template <typename T>
void pref_attach2(int npoints, int size, int dim, float radius, int iterations, mt19937 &mt, int time, bool mode) {
    int measurements = 10;
    vector<vector<int>> coord_list(npoints, vector<int>(Point<T>::Dim));
    vector<int> indices;
    vector<float> weights;
    vector<map<int, int>> mass_list(measurements);
    vector<map<int, int>>::iterator it = mass_list.begin();

    vector<vector<Neighbour<T, Point>>> neighbour_list(npoints);

    vector<Point<T>> points;
    points.reserve(npoints);

    Tree<int, Point> tree(points, float(size));
    gen_init_cond<T, Point>(points, coord_list, npoints, mt, size);
    tree.build(points);
    tree.find_neighbours(indices, neighbour_list, radius);
    tree.clear();
    for (int iteration = 0; iteration < iterations; iteration++) {
        for (int step = 0; step <= time; step++) {
            int rand_idx = static_cast<int>(mt() % static_cast<unsigned int>(points.size() - 1));
            while(points[rand_idx].mass == 0) rand_idx = static_cast<int>(mt() % static_cast<unsigned int>(points.size() - 1));
            Point<T>& point = points[rand_idx];
            for (int i = 0; i < neighbour_list[rand_idx].size(); i++){
                if (neighbour_list[rand_idx][i].point.mass != 0) {
                    weights.push_back(neighbour_list[rand_idx][i].point.mass/ pow(neighbour_list[rand_idx][i].distance, 2));
                }
            }
            discrete_distribution<int> dist(weights.begin(), weights.end());
            switch (mode) {
            case true:
                point.mass--;
                points[dist(mt)].mass ++;
                break;
            
            case false:
                points[dist(mt)].mass += point.mass;
                point.mass = 0;
                break;
            }
            weights.clear();
            cout << "\r" << string(30, ' ') << flush;
            cout << step/float(time) * 100 << "%" << flush;
            if (step == time/2){
                for (int i = 0; i < int(points.size()); i++) {
                    if (points[i].mass == 0) {
                        swap(points[i], points.back());
                        points.pop_back();
                        neighbour_list[i].clear();
                    }
                }
                tree.build(points);
                tree.find_neighbours(indices, neighbour_list, radius);
                tree.clear();
            }
            if (step != 0  && (step % (time/10)) == 0) {
                for (int i = 0; i < points.size(); i++) {
                    (*it)[points[i].mass] ++;
                }
                it++;
            }
        }
    }
    for (int i = 0; i < int(mass_list.size()); i++) {
        string name = "r" + to_string(int(radius)) + "_" + to_string(int(dim)) + "D_tree_" +to_string(i * time/10) +'_'+to_string(npoints/pow(size, dim)) + '_' + to_string(mode) +".txt";
        tree.output(name, mass_list[i], iterations, npoints, radius, size, dim, true);
    }
}

int main(){
    mt19937 mt(time(0));
    const unsigned long size = 50;
    const int npoints = 100000;
    const int time = 500000;
    const int dim = 3;
    float radius = 6;
    int iterations = 3;
    bool mode = true;
    

    if (dim != Point<int>::Dim) throw runtime_error("Dimension mismatch");
    //pref_growth<int>(npoints, size, dim, radius, iterations, mt);
    pref_attach<int>(npoints, size, dim, radius, iterations, mt , time, mode);
    return 1;
};