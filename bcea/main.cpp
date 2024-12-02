
#include "getArgs.hpp"
#include "BCE.h"
#include <string>
#include <chrono>
#include <iostream>

int main(int argc, char* argv[]) {
    std::string fPath;
    double mswbe_delta;
    unsigned int mswbe_left_thd;
    unsigned int mswbe_right_thd;
    unsigned int mswbe_reduction_strategy;
    unsigned int mswbe_sort_strategy;
   if(argc < 5){
        cout << "Too few arguments" << endl;
        return -1;
    }
    else{
        //for( int i = 0; i < argc; ++i )
        //    printf( "argv[%d]=%s\n", i, argv[i] );
        fPath = argv[1];
        mswbe_delta=atof(argv[2]);
        mswbe_left_thd=atoi(argv[3]);
        mswbe_right_thd=atoi(argv[4]);
        mswbe_reduction_strategy=atoi(argv[5]);
        mswbe_sort_strategy=atoi(argv[6]);
        
    }
   
    std::string order = "two";
    uint32_t ls = 1, rs = 1;
    

   
    std::cout << "file path: " << fPath << std::endl;
    std::cout << "graph order: " << order << std::endl;
    std::cout << "l " << ls << std::endl;
    std::cout << "r " << rs << std::endl;
    std::cout << "mswbe_delta："<<mswbe_delta<<std::endl;
    std::cout << "mswbe_left_thd："<<mswbe_left_thd<<std::endl;
    std::cout << "mswbe_right_thd："<<mswbe_right_thd<<std::endl;
    std::cout << "mswbe_reduction_strategy："<<mswbe_reduction_strategy<<std::endl;
    std::cout << "mswbe_sort_strategy："<<mswbe_sort_strategy<<std::endl;
    
    int flag=0;
    if(flag==1) {
        BCE bce(fPath, 0, order, ls, rs,mswbe_delta,mswbe_left_thd,mswbe_right_thd,mswbe_reduction_strategy,mswbe_sort_strategy);

        auto t1 = std::chrono::steady_clock::now();

        bce.run();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    }
    else {
        BCE bce(fPath, 1, order, ls, rs,mswbe_delta,mswbe_left_thd,mswbe_right_thd,mswbe_reduction_strategy,mswbe_sort_strategy);
        

        auto t1 = std::chrono::steady_clock::now();

        bce.run();
        removeDuplicates(bce.result_mswbe);
        bce.result_mswbe=reducePairs(bce.result_mswbe);
        //printResult(bce.result_mswbe);
        std::cout << "mswb number is:" << bce.result_mswbe.size()<<std::endl;
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        std::cout << "Memory usage : " << getMemoryUse() / (1024) << "(MB)" << endl;
    }


    return 0;
}

