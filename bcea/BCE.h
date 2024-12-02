#ifndef BCE_H
#define BCE_H

#include "biGraph.hpp"
#include "linearSetThree.hpp"
#include <string>
#include "mswbe.h"




class BCE {
private:
    biGraph * g;
    uint32_t ls, rs, lrs[2];
    std::string od;

    uint64_t totalMaximalCount = 0;

private:
    linearSetThree S[2];
    std::vector<std::vector<uint32_t>> ws;
    std::vector<uint32_t> realR[2];
    uint32_t okDegreeSide = 0;
    std::vector<uint32_t> deg;
    uint32_t realRSize[2];
    void bbranch(uint32_t);
    void bbranch2(uint32_t deep, uint32_t x[2], uint32_t c[2], uint32_t r[2]);

public:
unsigned int mswbe_left_thd,mswbe_right_thd,mswbe_reduction_strategy,mswbe_sort_strategy; 
    double mswbe_delta;
    BCE(const std::string & fPath, int mode, const std::string & order, uint32_t ls, uint32_t rs,unsigned int mswbe_delta,unsigned int mswbe_left_thd,unsigned int mswbe_right_thd,unsigned int mswbe_reduction_strategy,unsigned int mswbe_sort_strategy) {
        g = new biGraph(fPath, mode, order);
        
        printf("load graph\n");fflush(stdout);
        od = order;

        this->ls = ls;
        this->rs = rs;
        lrs[0] = ls;//?
        lrs[1] = rs;//?

        S[0].resize(g->n[0]);
        S[1].resize(g->n[1]);
        deg.resize(std::max(g->n[0], g->n[1]));//degree

        ws.resize(g->maxDu + g->maxDv);//
        realRSize[0] = realRSize[1] = 0;
        this->mswbe_delta=mswbe_delta;
        this->mswbe_left_thd=mswbe_left_thd;
        this->mswbe_right_thd=mswbe_right_thd;
        this->mswbe_reduction_strategy=mswbe_reduction_strategy;
        this->mswbe_sort_strategy=mswbe_sort_strategy;

    }
      
    std::vector< std::pair<std::vector<unsigned int>, std::vector<unsigned int> > > result_mswbe;
    void run();


};

#endif