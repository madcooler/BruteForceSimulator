//
//  main.cpp
//  
//
//  Created by Chi Wang on 04/12/15.
//  Copyright © 2015 Chi Wang. All rights reserved.
//


#include <iostream>
#include "BruteForceSimulator.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    
    BruteForceSimulator layer;
    
    //std::cout<<"Running..."<<"\n";
    //double start = omp_get_wtime( );
    layer.Calculate();
    //double end = omp_get_wtime( );
    //std::cout<<"Time cost £∫"<<end -start<<"\n";
    
    std::cout << "Hello, World!\n";
    return 0;
}
