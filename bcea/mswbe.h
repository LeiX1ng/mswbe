//
//  mswbe.h
//  
//
//  Created by 杨建业 on 2022/7/4.
//

#ifndef mswbe_h
#define mswbe_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <time.h>
#include <cassert>
#include <cmath>
#include <sys/resource.h>
#include <sys/time.h>
#include "naive.h"
#define COUNT_ONLY

//#define IS_DEBUGGING

//#define USE_SMALL_SIDE

#define INF 1000000

#define EPSILON 0.000001

#define mswbe_max(x,y) ((x)>(y)?(x):(y))
#define mswbe_min(x,y) ((x)<(y)?(x):(y))

using namespace std;

typedef unsigned long long lint;
typedef unsigned int ui;
typedef pair<ui, double> udp;
typedef pair<ui, ui> uup;
typedef pair<double, double> ddp;

class Tools{
public:
    static string integer_to_string(long long number);
};

class VertexDegree{
public:
    ui vertex;
    double degree;
    inline bool operator < (const VertexDegree &other) const {
        return degree < other.degree || (degree == other.degree && vertex < other.vertex);
    }

    VertexDegree();
    VertexDegree(ui v, double d);
    ~VertexDegree();
};

class Node{
public:
    vector<ui> L;
    vector<ui> R;
    vector<ui> CR;
    vector<ui> XR;
    double minW;
    double maxW;

    Node(vector<ui> L, vector<ui> R, vector<ui> CR, vector<ui> XR, double minW, double maxW);
    ~Node();
};

struct HashFunc {
    size_t operator() (const vector<ui>& key) const {
        std::hash<ui> hasher;
        size_t seed = 0;
        for (ui i : key) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct EqualFunc {
    bool operator() (vector<ui> const& a, vector<ui> const& b) const {
        for (ui i = 0; i < a.size(); i++) {
            if (a[i] != b[i])
                return 0;
        }
        return 1;
    }
};


class MSWBE{
public:

    MSWBE(string path, double delta_value, int left_size, int right_size, int reduction_strategy=1, int sort_strategy=1);
    ~MSWBE();
    void read_graph_txt(string file_path);
    void read_graph_txt(biclique bc);
    //void advanced_mswbe_old();
    void advanced_mswbe();
    void baseline_mswbe();
    void delete_memory();
    long getMemoryUse();

    vector< pair<vector<ui>, vector<ui> > > result_bicliques;
    lint result_count;

private:
    // the original graph
    ui n, n1, n2, m; //# of vertices, # of L-side vertices, # of R-side vertices, # of edges
    ui * pstart; //adjacent array
    udp * edges; //adjacent array
    int * degree; //vertex degree
    int * TMPdeg; //used in vertex reduction
    int * twoDeg; //used in vertex reduction
    int * side; //vertex side

    // the remaining graph
    ui r_n, r_n1, r_n2, r_m; //# of vertices, # of L-side vertices, # of R-side vertices, # of edges
    ui * r_pstart; //adjacent array
    udp * r_edges; //adjacent array
    int * r_degree; //vertex degree
    ui * oid; //original vertex id
    ui * nid; //new vertex id

    ui * tmp_r_pstart; //adjacent array
    udp * tmp_r_edges; //adjacent array
    int * tmp_r_degree; //vertex degree
    ui * tmp_oid; //original vertex id
    ui * tmp_nid; //new vertex id

    // the subgraph in enumeration procedure
    ui max_CS_deg; //used to construct Matrix
    double ** Matrix; //adjacent Matrix
    ui * trans; //mapping vertex id in Matrix
    bool * del_ver; //record vertex status
    int * inL; //in CR or not
    int * deg_inL; //degree in CR

    ddp * minmaxVec; // the min and max values of vertices in L

    // input parameters
    double delta; // weight similar threshold
    int left_thd, right_thd; // size constraints
    string data_file_path;
    int graph_reduction_strategy;
    int vertex_order_strategy;

    void core_based_reduction();

    void weight_based_reduction();

    int edgesSatisfied(vector<udp>& P, int size_thd);

    // get the two hop neighbors of v \in V among the neighbors of vertices in L
    vector<ui> get_two_hop_nei(ui v, vector<ui> &L, double lowW, double highW);

    // get the two hop neighbors of v \in V
    vector<ui> get_two_hop_nei(ui v);

    vector<ui> get_two_hop_nei(ui v, ui* &start_arr, udp* &edges_arr);

    vector<ui> get_two_hop_nei_with_range(ui v, ui* &start_arr, udp* &edges_arr, double &minW, double &maxW);

    void build_matrix_for_LCRXR(vector<ui> L, vector<ui> R, vector<ui> CR, vector<ui> XR);

    // minW and maxW denote the minimum and maximum edge weight in the biclique (clique_vertices_in_left, clique_vertices_in_right), respectively
    //void branch_and_bound(vector<ui> &L, vector<ui> &R, vector<ui> &CR, vector<ui> &XR, double minW, double maxW);

    void branch_and_bound(vector<ui> &L, vector<ui> &R, vector<ui> &CR, vector<ui> &XR);

    //
    //void level_based_search_old(Node node);

    void level_based_search(vector<udp> &P, ui vid);

    vector< pair<int, int> > find_sim_groups(vector<udp> &P, double d);

    vector< vector<ui> > collect_branches(vector< pair<ui, ddp> > &P, double d);

    void sort_vertices(vector<ui> &ids, ui strategy);


public:
    // tools
    void print(vector<ui> &vec);
    static void print_biclique(vector<ui> &L, vector<ui> &R);
    static void print_biclique(ofstream &fout, vector<ui> &L, vector<ui> &R);
    void print_biclique(vector<pair<vector<unsigned int>,vector<unsigned int>>>& callmswberesult,vector<ui> &L, vector<ui> &R,biclique& bc);
    void print_results(vector<pair<vector<unsigned int>,vector<unsigned int>>>& callmswberesult,biclique& bc);
    void print_results();
    void check_results();
    void print(vector<udp> &vec);
    void print(Node &node);
    void print(vector< pair<ui, ddp> > &P);
    static bool isSubsetOf(vector<ui> &Sub, vector<ui> &Sup);
    static bool isSameBiclique(pair<vector<ui>, vector<ui> > &b1, pair<vector<ui>, vector<ui> > &b2);
    static bool isSameVector(vector<ui> &vec1, vector<ui> &vec2);
    
    
};

void removeDuplicates(vector<pair<vector<unsigned int>, vector<unsigned int>>> &pairs);
bool contains(const vector<unsigned int> &a, const vector<unsigned int> &b);
vector<pair<vector<unsigned int>, vector<unsigned int>>> reducePairs(vector<pair<vector<unsigned int>, vector<unsigned int>>> pairs);

#endif /* mswbe_h */
