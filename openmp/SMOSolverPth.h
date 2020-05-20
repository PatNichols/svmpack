/*
 * SMOSolverPth.h
 *
 *  Created on: Jul 9, 2010
 *      Author: d3p708
 */

#ifndef SMOSOLVERPTH_H_
#define SMOSOLVERPTH_H_
#include "svm_utils.h"
#include "SVMOptions.h"
#include "SVMKernelMatrix.h"
#include "svm_traits.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <pthread.h>
#include "SVMKernelMatrixPth.h"
#include "CSMOThread.h"
#include "DSMOThread.h"
using namespace std;

namespace svmpack
{

namespace solv_run
{
template <class svm_real> static void * DSolverRun ( void *arg )
{
    DSMOThread<svm_real> *dsolver =
        reinterpret_cast< DSMOThread<svm_real> * > ( arg );
    dsolver->run();
    return 0x0;
}

template<class svm_real> static void * CSolverRun ( void *arg )
{
    CSMOThread<svm_real> *csolver =
        reinterpret_cast<CSMOThread<svm_real> *> ( arg );
    csolver->run();
    return 0x0;
}
}

template <class svm_real>
class SMOSolverPth
{
public:
    SMOSolverPth ( const SVMOptions<svm_real>& options_in ) :
            alpha ( 0x0 ), grad ( 0x0 ), dy ( options_in.getYAlphaPtr() ),  y ( 0x0 ), status ( 0x0 ),
            scal ( new svm_real[options_in.getNVectors() ] ),
            nvecs ( options_in.getNVectors() ),
            maxits ( options_in.getMaxIterations() ), fsum ( 0 ), bias ( 0 ),
            gap ( 0 ), eps ( options_in.getEpsilon() ),
            cost ( options_in.getCost() ), options ( options_in ) {
        try {
            alpha = new svm_real[nvecs];
            grad = new svm_real[nvecs];
            y = new int[nvecs];
            status = new int[nvecs];
        } catch ( exception& e ) {
            cerr << "could not allocate arrays for smo solver\n";
            exit ( EXIT_FAILURE );
        }
#ifdef _CHECK_DY_
        for ( size_t k = 0; k < nvecs; ++k ) {
            svm_real t=fabs(dy[k])-svm_real(1);
            if (t>eps) {
                cerr << " dy is out of whack ! " << dy[k] << " "<< k << endl;
            }
        }
#endif
        const svm_real zero ( 0 );
        for ( size_t k = 0; k < nvecs; ++k ) {
            alpha[k] = zero;
        }
        for ( size_t k = 0; k < nvecs; ++k ) {
            grad[k] = dy[k];
        }
        for ( size_t k = 0; k < nvecs; ++k ) {
            y[k] = ( dy[k] > zero ) ? 1 : -1;
        }
        for ( size_t k = 0; k < nvecs; ++k ) {
            status[k] = -1;
        }
    }
    ;

    ~SMOSolverPth() {
        delete[] status;
        delete[] y;
        delete[] grad;
        delete[] alpha;
    }
    ;

    svm_real getBias() const throw () {
        return bias;
    }
    ;
    svm_real getObjectiveFunction() const throw () {
        return fsum;
    }
    ;
    svm_real getGap() const throw () {
        return gap;
    }
    ;

    void findGap() throw () {
        const svm_real zero ( 0 );
        register svm_real asum ( 0 ), csum ( 0 );
        register int nsum ( 0 );
        fsum = 0;
        bias = 0;
        for ( size_t k = 0; k < nvecs; ++k ) {
            asum += alpha[k];
            fsum += alpha[k] * grad[k] * dy[k];
        }
        fsum = ( fsum + asum ) / 2;
        for ( size_t k = 0; k < nvecs; ++k ) {
            if ( status[k] )
                continue;
            bias += grad[k];
            nsum += 1;
        }
        if ( nsum )
            bias = -bias / nsum;
        for ( size_t k = 0; k < nvecs; ++k ) {
            csum += svmpack::fmax<svm_real> ( zero, ( ( grad[k] + bias ) * dy[k] ) );
        }
        csum *= cost;
        gap = ( asum + csum - fsum - fsum ) / ( asum + csum + svm_real ( 1 ) - fsum );
    }
    ;

    void outputModelFile() throw () {
        size_t nsv = 0;
        size_t nbnd = 0;
        for ( size_t k = 0; k < nvecs; ++k ) {
            if ( status[k] >= 0 ) {
                ++nsv;
                if ( status[k] > 0 )
                    ++nbnd;
            }
        }
        cerr << "# training vectors     = " << nvecs << endl;
        cerr << "# support vectors      = " << nsv << endl;
        cerr << "# bound support vectors= " << nbnd << endl;
        ofstream out ( options.getModelFileName().c_str() );
        if ( !out ) {
            cerr << "could not open file " << options.getModelFileName()
                 << " to output model \n";
            exit ( EXIT_FAILURE );
        }
        int itmp = static_cast<int> ( nsv );
        out.write ( ( char* ) &itmp, sizeof ( int ) );
        itmp = static_cast<int> ( options.getNFeatures() );
        out.write ( ( char* ) &itmp, sizeof ( int ) );
        itmp = static_cast<int> ( options.getKernelType() );
        out.write ( ( char* ) &itmp, sizeof ( int ) );
        itmp = static_cast<int> ( options.getKernelPower() );
        out.write ( ( char* ) &itmp, sizeof ( int ) );
        double dtmp = static_cast<double> ( options.getKernelCof1() );
        out.write ( ( char* ) &dtmp, sizeof ( double ) );
        dtmp = static_cast<double> ( options.getKernelCof2() );
        out.write ( ( char* ) &dtmp, sizeof ( double ) );
        dtmp = static_cast<double> ( bias );
        out.write ( ( char* ) &dtmp, sizeof ( double ) );
        itmp = ( options.scaleKernel() ) ? 1 : 0;
        out.write ( ( char* ) &itmp, sizeof ( int ) );
        for ( size_t k = 0; k < nvecs; ++k ) {
            if ( status[k] < 0 )
                continue;
            dtmp = dy[k] * scal[k] * alpha[k];
            out.write ( ( char* ) &dtmp, sizeof ( double ) );
        }
        if ( sizeof ( svm_real ) == sizeof ( double ) ) {
            size_t nfeat = options.getNFeatures();
            const svm_real *vp = options.getVectorsPtr();
            for ( size_t k = 0; k < nvecs; ++k ) {
                if ( status[k] < 0 ) continue;
                out.write ( ( char* ) ( vp + k * nfeat ), sizeof ( double ) *nfeat );
            }
        } else {
            size_t nfeat = options.getNFeatures();
            const svm_real *vp = options.getVectorsPtr();
            for ( size_t j = 0; j < nvecs; ++j ) {
                if ( status[j] < 0 ) continue;
                for ( size_t k = 0; k < nfeat; ++k ) {
                    dtmp = static_cast<double> ( vp[k+j*nfeat] );
                    out.write ( ( char* ) &dtmp, sizeof ( double ) );
                }
            }
        }
        out.close();
        ofstream out2 ( options.getOutputFileName().c_str() );
        if ( !out2 ) {
            cerr << "could not open file " << options.getOutputFileName()
                 << " to for output of scores \n";
            exit ( EXIT_FAILURE );
        }
        size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
        for ( size_t k = 0; k < nvecs; ++k ) {
            svm_real fx = ( dy[k] - grad[k] - bias );
            if ( fx > 0 ) {
                if ( y[k] > 0 )
                    ++ntp;
                else
                    ++nfp;
            } else {
                if ( y[k] > 0 )
                    ++ntn;
                else
                    ++nfn;
            }
            dtmp = static_cast<double> ( dy[k] );
            out2.write ( ( char* ) &dtmp, sizeof ( double ) );
            dtmp = static_cast<double> ( fx );
            out2.write ( ( char* ) &dtmp, sizeof ( double ) );
        }
        out2.close();
        analyze ( ntp, nfp, ntn, nfn );
    }
    ;


    void ctrain() {
        const size_t nth = options.getNThreads();
        pthread_t * pth = new pthread_t[nth];
        CSMOThread<svm_real> **cth =
            new CSMOThread<svm_real> * [nth];
        cth[0] = new CSMOThread<svm_real> ( options, cost, eps, alpha,
                                            grad, dy, y, status, nvecs, maxits, nth );
        memcpy ( scal, cth[0]->getScaleFactorPtr(), sizeof ( svm_real ) *nvecs );
        for ( size_t k = 1; k < nth; ++k ) {
            cth[k] = new CSMOThread<svm_real> ( options, **cth, k );
        }
        for ( size_t k = 0; k < nth; ++k ) {
            pthread_create ( ( pth + k ), 0x0, solv_run::CSolverRun<svm_real>, cth[k] );
        }
        for ( size_t k = 0; k < nth; ++k ) {
            pthread_join ( pth[k], 0x0 );
        }
        for ( size_t k = 0; k < nth; ++k ) delete cth[k];
        delete [] cth;
        delete [] pth;
    };

    void dtrain() {
        svm_stopwatch timer1;
        svm_stopwatch timer2;
        timer1.start();
        SVMKernelMatrix<svm_real> kmatrix ( options );
        memcpy ( scal, kmatrix.getScaleFactorPtr(), sizeof ( svm_real ) *nvecs );
        const svm_real *kmat = kmatrix.getKernelMatrixPtr();
        const size_t nth = options.getNThreads();
        pthread_t * pth = new pthread_t[nth];
        DSMOThread<svm_real> **dth =
            new DSMOThread<svm_real>*[nth];
        dth[0] = new DSMOThread<svm_real> ( cost, eps, alpha, grad,
                                            options.getYAlphaPtr(),
                                            kmat, y, status, nvecs, maxits, nth );
        for ( size_t k = 1; k < nth; ++k ) {
            dth[k] = new DSMOThread<svm_real> ( **dth, k );
        }
        timer2.start();
        for ( size_t k = 0; k < nth; ++k ) {
            pthread_create ( ( pth + k ), 0x0, solv_run::DSolverRun<svm_real>, dth[k] );
        }
        for ( size_t k = 0; k < nth; ++k ) {
            pthread_join ( pth[k], 0x0 );
        }
        timer2.stop();
        for ( size_t k = 0; k < nth; ++k ) delete dth[k];
        delete [] dth;
        delete [] pth;
        timer1.stop();
        cerr << "timer1 = " << timer1.elapsedTime() << " " << timer2.elapsedTime() << endl;
    };

    void train() throw() {
        svm_stopwatch timer;
        timer.start();
        size_t csize = options.getCacheSize();
        if ( csize > 1 ) ctrain();
        else dtrain();
        findGap();
        timer.stop();
        cerr << "training time = " << timer.elapsedTime() << endl;
    };
private:
    svm_real *alpha;
    svm_real *grad;
    const svm_real *dy;
    int *y;
    int *status;
    svm_real *scal;
    int step_return;
    size_t nvecs, maxits;
    svm_real fsum;
    svm_real bias;
    svm_real gap;
    svm_real eps, cost;
    const SVMOptions<svm_real>& options;
};

}

#endif /* SMOSOLVERPTH_H_ */
