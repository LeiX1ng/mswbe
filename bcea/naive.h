#ifndef NAIVE_H
#define NAIVE_H

#include <vector>
#include "biGraph.hpp"
#include <sys/resource.h>
//
// Created by Lei Xing on 2023/11/1
//
using namespace std;
struct weight_biclique{
    unsigned int u;
    unsigned int v;
    double weight;
};

class biclique
{
public:
    vector<unsigned int> vertices[2];
    unsigned int left_size;
    unsigned int right_size;
    vector<weight_biclique> wbc;
    unordered_map<unsigned int,unsigned int>leftMapindex;
    unordered_map<unsigned int,unsigned int>rightMapindex;
    unordered_map<unsigned int,unsigned int>indexMapleft;
    unordered_map<unsigned int,unsigned int>indexMapright;
    ~biclique();
    void get_wbc(biGraph& bg);
    void sorted_wbc();
    void print_wbc();
    void print_biclique();
};

void printResult(const std::vector<std::pair<std::vector<unsigned int>, std::vector<unsigned int>>>& result);
long getMemoryUse();

#endif