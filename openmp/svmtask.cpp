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
#include "svm_stopwatch.h"
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
    if ( options.getTask() == 0 ) {
        timer.clear();
        timer.start();
        if ( options.getNThreads() ) {
            svmpack::SMOSolverPth<SVMREAL> solver ( options );
            solver.train();
            solver.outputModelFile();
        } else {
#pragma omp parallel 
            {
                int nth = omp_get_num_threads();
#pragma omp single
                {
                    std::cerr << " # of openmp threads = " << nth << "\n"; 
                }                
            }
            svmpack::SMOSolver<SVMREAL> solver ( options );
            solver.train();
            solver.outputModelFile();
        }
        timer.stop();
        cerr << "time to train and write model = " << timer.elapsedTime() << "s \n";
    } else {
        timer.clear();
        timer.start();
        svmpack::SVMPredictor<SVMREAL> pred ( options );
        pred.predict();
        timer.stop();
        cerr << "time to predict = " << timer.elapsedTime() << "s \n";
    }
}
