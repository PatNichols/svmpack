//============================================================================
// Name        : 0.cpp
// Author      : Patrick Nichols
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "SVMOptions.h"
#include "svm_stopwatch.h"
#include "svm_utils.h"
using namespace std;

int main ( int argc, char **argv )
{
    typedef double SVMREAL;
    svmpack::svm_stopwatch timer;
    timer.start();
    svmpack::SVMOptions<SVMREAL> options ( argc, argv );
    timer.stop();
    cerr << "Time to get options " << timer.elapsedTime() << "s \n";
    cerr << options << endl;
    std::string filename = options.getDataFileName();
    timer.clear();
    timer.start();
    if (filename.find(".tdo")!=string::npos) {
        svmpack::tdo2libsvm(filename.c_str());
    }else{
        svmpack::libsvm2tdo(filename.c_str());
    }
    timer.stop();
    cerr << "translation time is " << timer.elapsedTime() << "\n";
    return EXIT_SUCCESS;
}
