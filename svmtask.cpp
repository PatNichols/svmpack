//============================================================================
// Name        : 0.cpp
// Author      : Patrick Nichols
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "SVMOptions.h"
#include "SMOSolver.h"
#include "SMOSolverPth.h"
#include "SVMPredictor.h"
using namespace std;

int main ( int argc, char **argv )
{
    typedef double SVMREAL;
    svmpack::SVMOptions<SVMREAL> options ( argc, argv );
    cerr << options << endl;
    if ( options.getTask() == 0 ) {
        if ( options.getNThreads() ) {
            svmpack::SMOSolverPth<SVMREAL> solver ( options );
            solver.train();
            solver.outputModelFile();
        } else {
            svmpack::SMOSolver<SVMREAL> solver ( options );
            solver.train();
            solver.outputModelFile();
        }
    } else {
        svmpack::SVMPredictor<SVMREAL> pred ( options );
        pred.predict();
    }
}
