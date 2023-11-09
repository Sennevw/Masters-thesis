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

template<class Point>
class Tree{
    public:
        Node* root_;
        vector<Point>* points_;
        float limit;

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
        
        void build(vector<Point>& points) {
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
            for (int i = 0; i < Point::Dim; i++){
                dist += pow( min(fabs(a.coords[i] - b.coords[i]), limit - fabs(a.coords[i] - b.coords[i])), 2);
            }
            return sqrt(dist);
        }

        bool same_coord(const Point& a, const Point& b){
            for (int i = 0; i < Point::Dim; i++){
                if (a.coords[i] != b.coords[i]) return false;
            }
            return true;
        }

        bool are_same(Node* trial, Point& center){
            if (same_coord((*points_)[trial->idx], center)) {
                (*points_)[trial->idx].mass += 1;
                center.mass = 0;
                (*points_).pop_back();
                return true;
            }
            return false;
        }

        bool insert(Point& point){
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

        bool insert_recurs(Node* node, Point& insert, int depth){
            const int axis = depth % Point::Dim;
            if (are_same(node, insert)){
                (*points_)[node->idx].mass += 1;
                return false;
            }
                
            if (insert.coords[axis] < (*points_)[node->idx].coords[axis]){
                if (node->left == nullptr){
                    node->left = new Node();
                    node->left->idx = points_->size() - 1;
                    node->left->axis = (depth + 1) % Point::Dim;
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
                    node->right->axis = (depth + 1) % Point::Dim;
                    return true;
                }
                else {
                    return insert_recurs(node->right, insert ,depth + 1);
                }
            }
            return true;
        }

        Node* search_parent(Point& point) {
            return search_parent_recurs(root_, point, 0);
        }

        Node* search_parent_recurs(Node* node, Point& point, int depth){
            if (node == nullptr) return nullptr;
            const int axis = depth % Point::Dim;
            if (are_same(node->left, point)) return node;
            else if (are_same(node->right, point)) return node;
            else if (point.coords[axis] < (*points_)[node->idx].coords[axis]) return search_parent_recurs(node->left, point, depth + 1);
            else return search_parent_recurs(node->right, point, depth + 1);
        }

        Node* search_replacement(Node* node){
            if (node->left == nullptr && node->right == nullptr) return nullptr;
            else if (node->right != nullptr) {
                node = node->right;
                while(node->left != nullptr) node = node->left;
                return node;
            }
            else if (node->left != nullptr) {
                node = node->left;
                while(node->right != nullptr) node = node->right;
                return node;
            }
            return nullptr;
        }

        void deletion(Point& point){
            Node* parent = search_parent(point);
            Node* node;
            if ((*points_)[node->idx].coords[parent->axis] < ((*points_)[parent->idx].coords[parent->axis])) node = parent->left;
            else node = parent->right;
            Node* replacement_ptr = search_replacement(node);
            Node replacement = *replacement_ptr;
            if (replacement_ptr->left != nullptr || replacement_ptr->right != nullptr) {
                deletion((*points_)[replacement_ptr->idx]);
            }
            else if (replacement_ptr->left == nullptr && replacement_ptr->right == nullptr) {
                replacement.left = node->left;
                replacement.right = node->right;
                replacement.axis = node->axis;
                if (are_same(parent->left, point)) parent->left = &replacement;
                else parent->right = &replacement;
                delete node;
            }
            return; 
        }

        bool radius_search(Point& center, double radius, vector<int>& indices, bool duplicates){
            bool flip = true;
            radius_search_recurs(root_, center, radius, indices, flip, duplicates);
            Point virtual_center = center;
            vector<vector<int>> coord_list;
            if (!flip) return flip;
            for (int i = 0; i < Point::Dim; i++){
                if (limit - center.coords[i] <= radius){
                    coord_list.push_back(center.coords);
                    coord_list.back()[i] = center.coords[i] - limit;
                    }
            }
            if (coord_list.size() > 1) {
                for (int i = 0; i < Point::Dim; i++) {
                    for (int j = 0; j < int(coord_list.size()); j++){
                        for (int k = 0; k < Point::Dim; k++){
                            if (coord_list[j][k] > 0 && limit - coord_list[j][k] <= radius){
                                coord_list.push_back(coord_list[j]);
                                coord_list.back()[k] = coord_list[j][k] - limit;
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
            else if (coord_list.size() == 1) {
                virtual_center.coords = coord_list[0];
                radius_search_recurs(root_, virtual_center, radius, indices, flip, false);
            }
            return flip;
        };
           

        void circ_radius_search_recurs(Node* trial, Point& center, double radius, vector<int>& indices, bool& flip, bool duplicates){
            if (trial == nullptr) return;
            if (!flip) return;
            if (duplicates && are_same(trial, center)) flip = false;
            if (!flip) return;
            const Point& trial_p = (*points_)[trial->idx];
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

        void radius_search_recurs(Node* trial, Point& center, double radius, vector<int>& indices, bool& flip, bool duplicates){
            if (trial == nullptr) return;
            if (!flip) return;
            if (duplicates && are_same(trial, center)) flip = false;
            if (!flip) return;
            const Point& trial_p = (*points_)[trial->idx];
            if (distance(trial_p, center) <= radius && !are_same(trial, center)){
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

        void create_random_int(vector<vector<int>> coord_list, int npoints, mt19937& mt){
            for (int i = 0; i < npoints; i++) {
                for (int j = 0; j < Point::Dim; j++){
                    coord_list[i][j] = static_cast<int>(mt() % static_cast<unsigned int>(limit));
                }
            }
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


// probeer dit eens met doubles
void pref_growth(int npoints, int size, int dim, float radius, int iterations, mt19937 &mt){
    vector<Point> points;
    points.reserve(npoints);
    vector<int> indices;
    vector<float> weights;
    Tree<Point> tree(points, float(size));
    map<int, int> mass_list;
    vector<vector<int>> coord_list(npoints, vector<int>(dim));

    for (int iteration = 0; iteration < iterations; iteration++) {
        tree.create_random_int(coord_list, npoints, mt);        
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
        for (int i = 0; i < npoints; i++) {
            points[i].mass = 0;
            points[i].coords = {};
        }
        points.clear();
        tree.points_ = &points;
    }
    cout << to_string(dim) + "D_tree_r"+to_string(radius) +'_'+to_string(npoints/pow(size, dim))+".txt" << endl;
    ofstream myfile;
    myfile.open(to_string(dim) + "D_tree_r" + to_string(radius) +'_'+to_string(npoints/pow(size, dim)) + ".txt");
    myfile << iterations << " " << npoints << " " << radius << endl;
    for (const auto& n: mass_list)
        myfile << n.first << " " << n.second << endl;
    myfile.close();
};

void pref_attach(int size, int npoints, int time, int dim, float radius, int iterations, mt19937 &mt){
    
    vector<vector<int>> coord_list(npoints, vector<int>(Point::Dim));

    vector<int> indices;
    vector<float> weights;
    map<int, int> mass_list;

    vector<Point> points;
    points.reserve(npoints);
    vector<Point*> current_points;
    current_points.reserve(npoints);

    Tree<Point> tree(points, float(size));
    

    
    for (int iteration = 0; iteration < iterations; iteration++) {
        tree.create_random_int(coord_list, npoints, mt);
        for (int i = 0; i < npoints; i++){
            points.push_back(Point(coord_list[i]));
            current_points.push_back(&points.back());
            if (!tree.insert(points.back())) {
                points.pop_back();
                current_points.pop_back();
            }
        }
        for (int steps = 0; steps <= time; steps++) {
            int rand_idx = mt() % points.size();
            Point& point = points[rand_idx];
            tree.radius_search(point, radius, indices, false);
            if (indices.size() > 0) {
                for (int j = 0; j < int(indices.size()); j++) {
                    weights.push_back(points[indices[j]].mass/ pow(tree.circ_distance(points[indices[j]], points.back()), 2));
                } 
                discrete_distribution<int> dist(weights.begin(), weights.end());
                points[indices[int(dist(mt))]].mass += 1;
                point.mass-= 1;
                if (point.mass == 0) {
                    tree.deletion(point);
                    if (rand_idx != current_points.size() - 1) {
                        swap(current_points[rand_idx], current_points.back());
                    }
                    current_points.pop_back();
                }
            }
            if (steps % 100 == 0) {
                cout << steps << endl;
                tree.iterate_tree(tree.root_, mass_list);
                cout << to_string(dim) + "D_tree_r"+to_string(radius) +'_'+to_string(npoints/pow(size, dim))+".txt" << endl;
                ofstream myfile;
                myfile.open(to_string(dim) + "D_tree_r" + to_string(radius) +'_'+to_string(npoints/pow(size, dim)) + ".txt");
                myfile << iterations << " " << npoints << " " << radius << endl;
                for (const auto& n: mass_list)
                    myfile << n.first << " " << n.second << endl;
                myfile.close();
                mass_list.clear();
                
            }
        }
        tree.clear();
        mass_list.clear();
        points.clear();
        current_points.clear();
        tree.points_ = &points;
        
    }
    
};

int main(){
    mt19937 mt(time(0));
    const unsigned long size = 50;
    const int npoints = 100000;
    const int time = 100;
    const int dim = 3;
    float radius = 3;
    int iterations = 500;

    if (dim != Point::Dim) throw runtime_error("Dimension mismatch");
    pref_growth(npoints, size, dim, radius, iterations, mt);
    return 1;
};




/*
int main(){
    mt19937 mt(time(0));
    const unsigned long size = 10;
    const int npoints = 13;
    const int dim = 2;
    vector<Point> points_2 = {  Point({1, 1}), 
                                Point({0, 0}), 
                                Point({1, 1}), 
                                Point({8, 9}), 
                                Point({9, 9}), 
                                Point({1, 9}), 
                                Point({3, 6}), 
                                Point({6, 2}), 
                                Point({9, 1}), 
                                Point({9, 2}),
                                Point({2, 0}),
                                Point({3, 6}),
                                Point({8, 9})};
 
    if (dim != Point::Dim) throw runtime_error("Dimension mismatch");
    vector<Point> points;
    points.reserve(npoints);
    vector<int> indices;
    vector<float> weights;
    
    Tree<Point> tree(points, int(size));
    
    for (int i = 0; i < npoints; i++) {
        points.push_back(points_2[i]);
        if (tree.radius_search(points.back(), 2, indices, true)) {
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
        indices.clear();
        
    }
   
    cout <<"what the fuck" << endl;
    for (vector<Point>::iterator it = points.begin(); it != points.end(); it++) {
        cout << it->coords[0] << " " << it->coords[1] << " " << it->coords[2] << " " << it->mass << endl;
    }
    return 1;
    
};
*/