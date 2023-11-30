#ifndef LINEARSETTHREE_HPP
#define LINEARSETTHREE_HPP

#include <cstdint>
#include <algorithm>
#include <cstdio>

class linearSetThree {
private:
    uint32_t * pSet = nullptr;


public:
    uint32_t x, c, r;
    uint32_t sz = 0;
    linearSetThree() {}
    linearSetThree(uint32_t size):sz(size), c(0), x(0), r(size) {
        vSet = new uint32_t[size];
        pSet = new uint32_t[size];

        for(uint32_t i = 0; i < size; i++) {
            vSet[i] = pSet[i] = i;
        }
    }
    void resize(uint32_t size) {
        sz = size;
        c = 0;
        x = 0;
        r = size;
        vSet = new uint32_t[size];
        pSet = new uint32_t[size];

        for(uint32_t i = 0; i < size; i++) {
            vSet[i] = pSet[i] = i;
        }
    }

    void reset() {
        c = 0;
        x = 0;
        r = sz;
    }

    ~linearSetThree() {
        delete [] vSet;
        delete [] pSet;
    }

    void swapByPos(uint32_t i, uint32_t j) {
        std::swap(vSet[i], vSet[j]);
        std::swap(pSet[vSet[i]], pSet[vSet[j]]);
    }

    uint32_t operator [] (uint32_t i) {
        // if(i >= sz) {
        //     printf("%u %u\n", i, sz);fflush(stdout);
        //     exit(-1);
        // }
        return  vSet[i];
    }

    uint32_t pos(uint32_t v) {
        return pSet[v];
    }

    bool CIsEmpty() { return c == r; }
    bool XIsEmpty() { return x == c; }

    uint32_t cSize() {
        return r - c;
    }

    uint32_t * begin() { return vSet; }

    void print() {
        if(x!=0){
            printf("yincang:");
            for(uint32_t i = 0; i < x; i++) {
                printf("%u ", vSet[i]);
            }
        }
        printf("X:");
        for(uint32_t i = x; i < c; i++) {
            printf("%u ", vSet[i]);
        }
        printf("C:");
        for(uint32_t i = c; i < r; i++) {
            printf("%u ", vSet[i]);
        }
        printf("R:");
        for(uint32_t i = r; i < sz; i++) {
            printf("%u ", vSet[i]);
        }
    }

    void print(uint32_t x, uint32_t c, uint32_t r) {
        printf("X:");
        for(uint32_t i = x; i < c; i++) {
            printf("%u ", vSet[i]);
        }
        printf("C:");
        for(uint32_t i = c; i < r; i++) {
            printf("%u ", vSet[i]);
        }
        printf("R:");
        for(uint32_t i = r; i < sz; i++) {
            printf("%u ", vSet[i]);
        }
    }

    uint32_t * vSet = nullptr;
};

#endif