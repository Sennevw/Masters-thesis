#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <array>

using namespace std;

class Point{
    public:
        vector<int> coords;
        int mass;
        static const int Dim = 2;

        Point(){};

        Point(vector<int> coord){
            this->coords = coord;
            this->mass = 1;
        }        
};

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
        vector<Point> points_;
        int limit;

        Tree(): root_(nullptr){};

        Tree(vector<Point>& points, int lim): root_(nullptr), limit(lim){build(points);}

        
        ~Tree() { clear(); }
        
        void build(const vector<Point> &points) {
            clear();
            if (points[0].mass == 0) return;
            vector<int> indices(points.size());
            iota(indices.begin(), indices.end(), 0);
            root_ = build_recurs(indices.data(), (int)points.size(), 0);
        }

        void clear() {
            clear_recurs(root_);
            root_ = nullptr;
            points_.clear();
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
                return points_[lhs].coords[axis] < points_[rhs].coords[axis];
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
                dist += pow( min(abs(a.coords[i] - b.coords[i]), limit - abs(a.coords[i] - b.coords[i])), 2);
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
                root_->idx = points_.size();
                root_->axis = 0 % Point::Dim;
                points_.push_back(point);
                return true;
            }
            else {
                return insert_recurs(root_, point, 0);
            }
        }

        bool insert_recurs(Node* node, const Point& insert, int depth){
            const int axis = depth % Point::Dim;
            if (insert.coords[axis] < points_[node->idx].coords[axis]){
                if (node->left == nullptr){
                    node->left = new Node();
                    node->left->idx = points_.size();
                    node->left->axis = depth % Point::Dim;
                    points_.push_back(insert);
                    return true;
                }
                else {
                    return insert_recurs(node->left, insert, depth + 1);
                }
            }
            else {
                if (node->right == nullptr){
                    node->right = new Node();
                    node->right->idx = points_.size();
                    node->right->axis = depth % Point::Dim;
                    points_.push_back(insert);
                    return true;
                }
                else {
                    return insert_recurs(node->right, insert ,depth + 1);
                }
            } 
        }

        bool radius_search(const Point& center, double radius, vector<int>& indices, bool duplicates = false){
            bool flip = true;
            radius_search_recurs(root_, center, radius, indices, duplicates, bool &flip);
            return flip;
        }

        void radius_search_recurs(Node* trial, const Point& center, double radius, vector<int>& indices, bool duplicates = false, bool& flip){
            if (trial == nullptr) return 0;
            if (duplicates) {
                if (are_same(center, points_[trial->idx])) {
                    points_[trial->idx].mass += 1;
                    points_.push_back(center);
                    flip = false;
                    return;
                }
            }
            const Point& trial_p = points_[trial->idx];
            if (!duplicates) {
                if (circ_distance(trial_p, center) <= radius && !are_same(trial_p, center)){
                    indices.push_back(trial->idx);
                }
            }
            const int dir = center.coords[trial->axis] < trial_p.coords[trial->axis] ? 0 : 1;
            if (!dir) {
                radius_search_recurs(trial->left, center, radius, indices);
            }
            else {
                radius_search_recurs(trial->right, center, radius, indices);
            }
            
            const int axis = trial->axis;
            const double diff = fabs(center.coords[axis] - trial_p.coords[axis]);
            if (diff < radius) {
                radius_search_recurs(dir ? trial->left : trial->right, center, radius, indices);
            }
            return;
        }
};

void test_full_list(Tree<Point> &tree, vector<Point> &points, int npoints, int size, int dim, mt19937& mt) {
    for (int i = 0; i < npoints; i++) {
        vector<int> coords(dim);
        for (int j = 0; j < dim; j++){
            coords[j] = mt() % size;
        }
        points[i] = Point(coords);
        
        tree.insert(points[i]);
        cout << points[i].coords[0] << " " << points[i].coords[1] << endl;
    }
    cout << endl;
    vector<int> indices;
    tree.radius_search(points[1], 4, indices);
    cout << "printing indices" << endl;
    for (int i = 0; i < indices.size(); i++){
        cout << points[indices[i]].coords[0] << " " << points[indices[i]].coords[1] << endl;
    }
}

void test_gen_step(Tree<Point> &tree, vector<Point> &points, int npoints, int size, int dim, mt19937& mt) {
    vector<Point> points_2 = {  Point({0, 0}), 
                                Point({0, 1}), 
                                Point({1, 0}), 
                                Point({1, 1}), 
                                Point({0, 0}), 
                                Point({0, 1}), 
                                Point({1, 0}), 
                                Point({1, 1}), 
                                Point({0, 0}), 
                                Point({0, 1})};
    for (int i = 0; i < npoints; i++) {

    }
}

int main(){
    mt19937 mt(time(0));
    const unsigned long size = 10;
    const int npoints = 10;
    const int dim = 2;

    if (dim != Point::Dim) throw runtime_error("Dimension mismatch");

    
    vector<Point> points(npoints);
    Tree<Point> tree(points, int(size));


}