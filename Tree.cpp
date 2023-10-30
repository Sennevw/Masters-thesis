#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <array>
#include <iostream>
using namespace std;

struct Node {
    Node* left;
    Node* right;
    int idx;
    int axis;
    Node(): idx(-1), axis(-1){this->left = this->right =nullptr;};
};

template<class Point>
class Tree{
    public:
        Node* root_;
        vector<Point>* points_;
        int limit;

        Tree(int lim): root_(nullptr), points_(nullptr), limit(lim){};

        Tree(vector<Point>& points, int lim): root_(nullptr), limit(lim){
            if (!points.empty()) {
                points_ = &points;
                build(points);
            }
            else {
                points_ = &points;
            }
        }

        ~Tree() { clear(); }
        
        void build(const vector<Point> &points) {
            clear();
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
            const int axis = depth % Point::Dim;
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

        static double distance(const Point& a, const Point& b){
            double dist = 0;
            for (int i = 0; i < Point::Dim; ++i){
                dist += (a.coords[i] - b.coords[i]) * (a.coords[i] - b.coords[i]);
            }
            return sqrt(dist);
        }

        static double pow_distance(const Point& a, const Point& b){
            double dist = 0;
            for (int i = 0; i < Point::Dim; ++i){
                dist += (a.coords[i] - b.coords[i]) * (a.coords[i] - b.coords[i]);
            }
            return dist;
        }

        double circ_distance(const Point& a, const Point& b) {
            double dist = 0;
            for (int i = 0; i < Point::Dim; ++i){
                dist += pow( min(fabs(a.coords[i] - b.coords[i]), limit - fabs(a.coords[i] - b.coords[i])), 2);
            }
            return sqrt(dist);
        }

        bool are_same(const Point& a, const Point& b){
            for (int i = 0; i < Point::Dim; i++){
                if (a.coords[i] != b.coords[i]) return false;
            }
            return true;
        }

        bool insert(const Point& point){
            if (root_ == nullptr) {
                root_ = new Node();
                root_->idx = 0;
                root_->axis = 0 % Point::Dim;
                return true;
            }
            else {
                return insert_recurs(root_, point, 0);
            }
        }

        bool insert_recurs(Node* node, const Point& insert, int depth){
            const int axis = depth % Point::Dim;
            if (insert.coords[axis] < (*points_)[node->idx].coords[axis]){
                if (node->left == nullptr){
                    node->left = new Node();
                    node->left->idx = points_->size() - 1;
                    node->left->axis = depth % Point::Dim;
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
                    node->right->axis = depth % Point::Dim;
                    return true;
                }
                else {
                    return insert_recurs(node->right, insert ,depth + 1);
                }
            }
            return true;
        }

        bool radius_search(Point& center, double radius, vector<int>& indices, bool duplicates = false){
            bool flip = true;
            radius_search_recurs(root_, center, radius, indices, duplicates, flip);
            return flip;
        }

        void radius_search_recurs(Node* trial, Point& center, double radius, vector<int>& indices, bool& flip, bool duplicates = false){
            if (trial == nullptr) return;
            if (duplicates) {
                if (are_same(center, (*points_)[trial->idx])) {
                    (*points_)[trial->idx].mass += 1;
                    center.mass = 0;
                    (*points_).pop_back();
                    indices.clear();
                    flip = false;
                    return;
                }
            }
            const Point& trial_p = (*points_)[trial->idx];
            if (circ_distance(trial_p, center) <= radius && !are_same(trial_p, center)){
                indices.push_back(trial->idx);
            }
            const int dir = center.coords[trial->axis] < trial_p.coords[trial->axis] ? 0 : 1;
            if (!dir) {
                radius_search_recurs(trial->left, center, radius, indices, duplicates);
            }
            else {
                radius_search_recurs(trial->right, center, radius, indices, duplicates);
            }
            
            const int axis = trial->axis;
            const double diff = fabs(center.coords[axis] - trial_p.coords[axis]);
            if (diff < radius || ((limit - diff) < radius)) {
                radius_search_recurs(dir ? trial->left : trial->right, center, radius, indices, duplicates);
            }
            return;
        }

        void iterate_tree(Node* node, vector<int>& mass_list){
            if (node == nullptr) return;
            mass_list[(*points_)[node->idx].mass] += 1;
            iterate_tree(node->left, mass_list);
            iterate_tree(node->right, mass_list);
        }
};

class Point{
    public:
        vector<int> coords;
        int mass;
        static const int Dim = 3;

        Point(){};

        Point(vector<int> coord){
            this->coords = coord;
            this->mass = 1;
        }        
};


int main(){
    mt19937 mt(time(0));
    const unsigned long size = 50;
    const int npoints = 100000;
    const int dim = 3;
    float radius = 5;
    int iterations = 500;
    vector<int> mass_list(npoints);


    if (dim != Point::Dim) throw runtime_error("Dimension mismatch");
    vector<Point> points;
    points.reserve(npoints);
    vector<int> indices;
    vector<float> weights;
    Tree<Point> tree(points, int(size));

    vector<vector<int>> coord_list(npoints);

    for (int iteration = 0; iteration < iterations; iteration++) {
        vector<vector<int>> coord_list(npoints);
        for (int i = 0; i < npoints; i++) {
            vector<int> coord(dim);
            for (int j = 0; j < dim; j++){
                coord[j] = mt() % int(size);
            }
            coord_list[i] = coord;
        }
        cout << iteration << endl;
        for (int i = 0; i < npoints; i++) {
            points.push_back(Point(coord_list[i]));
            if (tree.radius_search(points.back(), radius, indices, true)) {
                if (indices.size() > 0) {
                    for (int j = 0; j < int(indices.size()); j++) {
                        weights.push_back(points[indices[j]].mass/ pow(tree.circ_distance(points[indices[j]], points.back()), 2));
                    } 
                    discrete_distribution<int> dist(weights.begin(), weights.end());
                    points[indices[int(dist(mt))]].mass += 1;
                    points.pop_back();
                    weights.clear();
                    
                }
                else {
                    tree.insert(points.back());
                }
            }
            indices.clear();
        }
        tree.iterate_tree(tree.root_, mass_list);
        tree.clear();
        points.clear();
        coord_list.clear();
        tree.points_ = &points;
    }
    ofstream myfile;
    myfile.open("tree_r" + to_string(radius) +'_'+to_string(npoints/pow(size, 3)) + ".txt");
    myfile << iterations << " " << npoints << " " << radius << endl;
    for (int i = 0; i < npoints; i++){
        myfile << i << " " << mass_list[i] << endl;
    }
    return 1;
};
 /*
int main(){
    mt19937 mt(time(0));
    const unsigned long size = 10;
    const int npoints = 12;
    const int dim = 2;
    vector<Point> points_2 = {  Point({1, 1}), 
                                Point({0, 0}), 
                                Point({0, 1}), 
                                Point({8, 9}), 
                                Point({9, 9}), 
                                Point({1, 9}), 
                                Point({3, 6}), 
                                Point({6, 2}), 
                                Point({9, 1}), 
                                Point({9, 2}),
                                Point({2, 0}), 
                                Point({3, 6})};
 
    if (dim != Point::Dim) throw runtime_error("Dimension mismatch");
    vector<Point> points;
    points.reserve(npoints);
    vector<int> indices;
    vector<float> weights;
    
    Tree<Point> tree(points, int(size));
    
    for (int i = 0; i < npoints; i++) {
        points.push_back(points_2[i]);
        tree.radius_search(points.back(), 3, indices, true);
        if (indices.size() > 0) {
            for (int j = 0; j < indices.size(); j++) {
                weights.push_back(points[indices[j]].mass/ tree.pow_distance(points[indices[j]], points.back()));
            } 
            discrete_distribution<int> dist(weights.begin(), weights.end());
            int random = indices[dist(mt)];
            points[random].mass += 1;
            points.pop_back();
            weights.clear();
            indices.clear();
        }
        else {
            tree.insert(points.back());
        }
        
    }
   
    cout <<"what the fuck" << endl;
    for (vector<Point>::iterator it = points.begin(); it != points.end(); it++) {
        cout << it->coords[0] << " " << it->coords[1] << " " << it->mass << endl;
    }
    return 1;
};
*/