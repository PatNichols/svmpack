/*
 * SVMPredictor.h
 *
 *  Created on: Jul 8, 2010
 *      Author: d3p708
 */
#ifndef SVMPREDICTOR_H_
#define SVMPREDICTOR_H_
#include "SVMOptions.h"
#include "SVMKernelEvaluator.h"
#include "svm_utils.h"
#include "svm_traits.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <pthread.h>
using namespace std;

namespace svmpack
{

template<class svm_real>
class SVMPredictor
{
public:
    SVMPredictor ( const SVMOptions<svm_real>& options_in ) :
            options ( options_in ) {
    }
    ;

    ~SVMPredictor() {
    }
    ;

    void predict() throw ();

    void predict_pth() throw ();

private:
    const SVMOptions<svm_real>& options;
};


template<class svm_real>
class SVMProductEvaluator
{
public:
    SVMProductEvaluator ( const SVMOptions<svm_real>& options ) :
            vecs ( options.getVectorsPtr() ),
            scal ( options.getYAlphaPtr() ),
            c1 ( options.getKernelCof1() ),
            c2 ( options.getKernelCof2() ),
            nvecs ( options.getNVectors() ),
            nfeat ( options.getNFeatures() ),
            ktype ( options.getKernelType() ),
            kpow ( options.getKernelPower() ) {
        if ( nvecs == 0 || nfeat == 0 || vecs == 0 ) {
            cerr << "error in input to SVMKernelEvaluator!\n";
            cerr << "nvectors  = " << nvecs << endl;
            cerr << "nfeatures  = " << nfeat << endl;
            cerr << "vector ptr= " << vecs << endl;
            exit ( EXIT_FAILURE );
        }
        if ( ktype < 0 || ktype > 3 ) {
            cerr << " unknown kernel type passed to SVMKernelEvaluator!\n";
            exit ( EXIT_FAILURE );
        }
        if ( ktype == 2 ) {
            c1 = -c1;
        }
    }
    ;

    ~SVMProductEvaluator() {
    };

    inline const svm_real *getScaleFactorPtr() const throw() {
        return scal;
    };

    inline svm_real getScaleFactor ( size_t ix ) const throw() {
        return scal[ix];
    };

    inline svm_real eval ( const svm_real * __restrict__ v1,
                           const svm_real * __restrict__ v2 ) const throw () {
        switch ( ktype ) {
        case 0:
            return eval0 ( v1, v2 );
        case 1:
            return eval1 ( v1, v2 );
        case 2:
            return eval2 ( v1, v2 );
        case 3:
            return eval3 ( v1, v2 );
        }
        return eval2 ( v1, v2 );
    }
    ;

    inline svm_real eval ( const svm_real * __restrict__ v1 ) const throw () {
        switch ( ktype ) {
        case 0:
            return eval0 ( v1, v1 );
        case 1:
            return eval1 ( v1, v1 );
        case 2:
            return svm_real ( 1 );
        case 3:
            return eval3 ( v1, v1 );
        }
        return svm_real ( 1 );
    }
    ;

    inline svm_real genSum ( const svm_real * __restrict__ v1 ) const throw () {
        register svm_real sum ( 0 );
        const svm_real * v2 = vecs;
        switch ( ktype ) {
        case 0:
            for ( size_t k = 0; k < nvecs; ++k, v2 += nfeat ) {
                sum += eval0 ( v1, v2 ) * scal[k];
            }
            return sum;
        case 1:
            for ( size_t k = 0; k < nvecs; ++k, v2 += nfeat ) {
                sum += eval1 ( v1, v2 ) * scal[k];
            }
            return sum;
        case 2:
            for ( size_t k = 0; k < nvecs; ++k, v2 += nfeat ) {
                sum += eval2 ( v1, v2 ) * scal[k];
            }
            return sum;
        case 3:
            for ( size_t k = 0; k < nvecs; ++k, v2 += nfeat ) {
                sum += eval3 ( v1, v2 ) * scal[k];
            }
            return sum;
        }
        return 0;
    }
    ;
private:
    const svm_real *vecs;
    const svm_real *scal;
    svm_real c1, c2;
    size_t nvecs, nfeat, ktype, kpow;


    inline svm_real eval0 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
        register svm_real s ( 0 );
        for ( register size_t k = 0; k < nfeat; ++k )
            s += v1[k] * v2[k];
        return s;
    }
    ;

    inline svm_real eval1 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
        register svm_real s ( 0 );
        for ( register size_t k = 0; k < nfeat; ++k )
            s += v1[k] * v2[k];
        return svmpack::powi<svm_real> ( ( s * c1 + c2 ), kpow );
    }
    ;

    inline svm_real eval3 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
        register svm_real s ( 0 );
        for ( register size_t k = 0; k < nfeat; ++k )
            s += v1[k] * v2[k];
        return tanh ( ( c1 * s + c2 ) );
    }
    ;

    inline svm_real eval2 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
        register svm_real t, s ( 0 );
        for ( register size_t k = 0; k < nfeat; ++k ) {
            t = v1[k] - v2[k];
            s += t * t;
        }
        s *= c1;
        return exp(s);
    }
    ;

};


template<class svm_real> struct predictor_thread {
    const SVMOptions<svm_real>& options;
    const SVMProductEvaluator<svm_real> kfun;
    size_t m_nvecs;
    size_t *arr;
    size_t tid, nth;

    predictor_thread ( const SVMOptions<svm_real>& options_in,
                       size_t thread_id_in, size_t *arr_in ) :
            options ( options_in ), kfun ( options_in ), m_nvecs ( options.getNVectors() ),
            arr ( arr_in ), tid ( thread_id_in ), nth ( options.getNThreads() ) {
    }
    ;

    ~predictor_thread() {
    }
    ;

    inline void predict() {
        ifstream in ( options.getDataFileName().c_str() );
        if ( !in ) {
            cerr << "could not open data file " << options.getDataFileName()
                 << "with data to predict\n";
            exit ( EXIT_FAILURE );
        }
        int itmp;
        in.read ( ( char* ) &itmp, sizeof ( int ) );
        const size_t d_nvecs = static_cast<size_t> ( itmp );
        in.read ( ( char* ) &itmp, sizeof ( int ) );
        const size_t d_nfeat = static_cast<size_t> ( itmp );
        size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
        const svm_real bias = options.getBias();
        const svm_real tau = svm_traits<svm_real>::tau();
        const bool scale = options.scaleKernel();
        double *dlabs;
        double *dvec;
        try {
            dlabs = new double[d_nvecs];
            dvec = new double[d_nfeat];
        } catch ( exception& e ) {
            cerr << " could not allocate array in " << __FILE__ << "("
                 << __LINE__ << ")\n";
            exit ( EXIT_FAILURE );
        }
        in.read ( ( char* ) dlabs, sizeof ( double ) * d_nvecs );
        size_t bsz = d_nvecs / nth;
        size_t xsz = d_nvecs % nth;
        size_t off1, off2;
        if ( tid < xsz ) {
            off1 = bsz * tid + tid;
            off2 = bsz + 1 + off1;
        } else {
            off1 = bsz * tid + xsz;
            off2 = bsz + off1;
        }
        size_t file_off1 = sizeof ( double ) * off1 * d_nfeat;
        in.seekg ( file_off1, ios::cur );
        svm_real * labs;
        svm_real * vec;
        try {
            labs = new svm_real[d_nvecs];
            vec = new svm_real[d_nfeat];
        } catch ( exception& e ) {
            cerr << " could not allocate array in " << __FILE__ << "("
                 << __LINE__ << ")\n";
            exit ( EXIT_FAILURE );
        }
        for ( size_t k = 0; k < d_nvecs; ++k ) {
            labs[k] = static_cast<svm_real> ( dlabs[k] );
        }
        for ( size_t k = off1; k < off2; ++k ) {
            in.read ( ( char* ) dvec, sizeof ( double ) * d_nfeat );
            for ( size_t l = 0; l < d_nfeat; ++l )
                vec[l] = static_cast<svm_real> ( dvec[l] );
            svm_real fx = kfun.genSum ( vec );
            if ( scale ) {
                svm_real sx = kfun.eval ( vec );
                if ( sx > tau )
                    sx = 1 / sqrt ( sx );
                else
                    sx = 1;
                fx *= sx;
            }
            fx -= bias;
            if ( fx > 0 ) {
                if ( labs[k] > 0 )
                    ++ntp;
                else
                    ++nfp;
            } else {
                if ( labs[k] > 0 )
                    ++ntn;
                else
                    ++nfn;
            }
        }
        in.close();
        arr[0] = ntp;
        arr[1] = nfp;
        arr[2] = ntn;
        arr[3] = nfn;
    }
    ;
};

template <> inline void predictor_thread<double>::predict()
{
    ifstream in ( options.getDataFileName().c_str() );
    if ( !in ) {
        cerr << "could not open data file " << options.getDataFileName()
             << "with data to predict\n";
        exit ( EXIT_FAILURE );
    }
    int itmp;
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    const size_t d_nvecs = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    const size_t d_nfeat = static_cast<size_t> ( itmp );
    size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
    const double bias = options.getBias();
    const double tau = svm_traits<double>::tau();
    const bool scale = options.scaleKernel();
    double *dlabs;
    double *dvec;
    try {
        dlabs = new double[d_nvecs];
        dvec = new double[d_nfeat];
    } catch ( exception& e ) {
        cerr << " could not allocate array in " << __FILE__ << "("
             << __LINE__ << ")\n";
        exit ( EXIT_FAILURE );
    }
    in.read ( ( char* ) dlabs, sizeof ( double ) * d_nvecs );
    size_t bsz = d_nvecs / nth;
    size_t xsz = d_nvecs % nth;
    size_t off1, off2;
    if ( tid < xsz ) {
        off1 = bsz * tid + tid;
        off2 = bsz + 1 + off1;
    } else {
        off1 = bsz * tid + xsz;
        off2 = bsz + off1;
    }
    size_t file_off1 = sizeof ( double ) * off1 * d_nfeat;
    in.seekg ( file_off1, ios::cur );
    for ( size_t k = off1; k < off2; ++k ) {
        in.read ( ( char* ) dvec, sizeof ( double ) * d_nfeat );
        double fx = 0;
        fx = kfun.genSum ( dvec );
        if ( scale ) {
            double sx = kfun.eval ( dvec );
            if ( sx > tau )
                sx = 1. / sqrt ( sx );
            else
                sx = 1.;
            fx *= sx;
        }
        fx -= bias;
        if ( fx > 0 ) {
            if ( dlabs[k] > 0 )
                ++ntp;
            else
                ++nfp;
        } else {
            if ( dlabs[k] > 0 )
                ++ntn;
            else
                ++nfn;
        }
    }
    arr[0] = ntp;
    arr[1] = nfp;
    arr[2] = ntn;
    arr[3] = nfn;
};


template<class svm_real>
static void * predictor_trun ( void * args )
{
    predictor_thread<svm_real> *t = reinterpret_cast< predictor_thread<svm_real>* > ( args );
    t->predict();
    return 0x0;
}

template<class svm_real>
inline void svmpack::SVMPredictor<svm_real>::predict_pth() throw ()
{
    size_t nth = options.getNThreads();
    size_t *arr = new size_t[nth * 4];
    size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
    predictor_thread<svm_real> **t = new predictor_thread<svm_real>*[nth];
    for ( size_t k = 0; k < nth; ++k ) {
        t[k] = new predictor_thread<svm_real> ( options, k, (arr + 4 * k) );
    }
    pthread_t *pth = new pthread_t[nth];
    for ( size_t k = 0; k < nth; ++k ) {
        pthread_create ( pth + k, 0x0, predictor_trun<svm_real> , t[k] );
    }
    void *ret_val;
    for ( size_t k = 0; k < nth; ++k ) {
        int err = pthread_join ( pth[k], &ret_val );
        if ( err ) {
            cerr << " thread " << k << "return error " << err << endl;
        }
    }
    for ( size_t k = 0; k < nth; ++k ) {
        ntp += arr[k * 4];
        nfp += arr[k * 4 + 1];
        ntn += arr[k * 4 + 2];
        nfn += arr[k * 4 + 3];
    }
    analyze ( ntp, nfp, ntn, nfn );
    delete[] pth;
    for ( size_t k = 0; k < nth; ++k )
        delete t[k];
    delete[] t;
    delete[] arr;
}

template<class svm_real> inline
void svmpack::SVMPredictor<svm_real>::predict() throw ()
{
    svm_stopwatch timer;
    timer.clear();
    timer.start();
    size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
    if ( options.getNThreads() != 0 ) {
        this->predict_pth();
        timer.stop();
        cerr << "prediction time  = " << timer.elapsedTime() << endl;
        return;
    }
    SVMProductEvaluator<svm_real> kfun ( options );
    ifstream in ( options.getDataFileName().c_str() );
    if ( !in ) {
        cerr << "could not open data file " << options.getDataFileName()
             << "with data to predict\n";
        exit ( EXIT_FAILURE );
    }
    int itmp;
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    size_t d_nvecs = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    size_t d_nfeat = static_cast<size_t> ( itmp );
    const svm_real bias = options.getBias();
    const bool scale = options.scaleKernel();
    const svm_real tau = svm_traits<svm_real>::tau();
    double *dlabs;
    double *dvec;
    svm_real *labs;
    svm_real *vec;
    try {
        dlabs = new double[d_nvecs];
        labs = new svm_real[d_nvecs];
        dvec = new double[d_nfeat];
        vec = new svm_real[d_nfeat];
    } catch ( exception& e ) {
        cerr << " could not allocate array in " << __FILE__ << "("
             << __LINE__ << ")\n";
        exit ( EXIT_FAILURE );
    }
    in.read ( ( char* ) dlabs, sizeof ( double ) * d_nvecs );
    for ( size_t k = 0; k < d_nvecs; ++k ) {
        labs[k] = static_cast<svm_real> ( dlabs[k] );
    }
    for ( size_t k = 0; k < d_nvecs; ++k ) {
        in.read ( ( char* ) dvec, sizeof ( double ) * d_nfeat );
        for ( size_t l = 0; l < d_nfeat; ++l )
            vec[l] = static_cast<svm_real> ( dvec[l] );
        svm_real fx = kfun.genSum ( vec );
        if ( scale ) {
            svm_real sx = kfun.eval ( vec );
            if ( sx > tau )
                sx = 1 / sqrt ( sx );
            else
                sx = 1;
            fx *= sx;
        }
        fx -= bias;
        if ( fx > 0 ) {
            if ( labs[k] > 0 )
                ++ntp;
            else
                ++nfp;
        } else {
            if ( labs[k] > 0 )
                ++ntn;
            else
                ++nfn;
        }
    }
    in.close();
    analyze ( ntp, nfp, ntn, nfn );
        timer.stop();
        cerr << "prediction time  = " << timer.elapsedTime() << endl;
        return;
}

template<> inline
void svmpack::SVMPredictor<double>::predict() throw ()
{
    svm_stopwatch timer;
    
    timer.start();
    if ( options.getNThreads() != 0 ) {
        this->predict_pth();
        timer.stop();
        cerr << "prediction time  = " << timer.elapsedTime() << endl;
        return;
    }
    ifstream in ( options.getDataFileName().c_str() );
    if ( !in ) {
        cerr << "could not open data file " << options.getDataFileName()
             << "with data to predict\n";
        exit ( EXIT_FAILURE );
    }
    int itmp;
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    size_t d_nvecs = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    size_t d_nfeat = static_cast<size_t> ( itmp );
    size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
    SVMProductEvaluator<double> kfun ( options );
    const double bias = options.getBias();
    const bool scale = options.scaleKernel();
    const double tau = svm_traits<double>::tau();
    double *dlabs;
    double *dvec;
    try {
        dlabs = new double[d_nvecs];
        dvec = new double[d_nfeat];
    } catch ( exception& e ) {
        cerr << " could not allocate array in " << __FILE__ << "("
             << __LINE__ << ")\n";
        exit ( EXIT_FAILURE );
    }
    in.read ( ( char* ) dlabs, sizeof ( double ) * d_nvecs );
    for ( size_t k = 0; k < d_nvecs; ++k ) {
        in.read ( ( char* ) dvec, sizeof ( double ) * d_nfeat );
        double fx = kfun.genSum ( dvec );
        if ( scale ) {
            double sx = kfun.eval ( dvec );
            if ( sx > tau )
                sx = 1 / sqrt ( sx );
            else
                sx = 1;
            fx *= sx;
        }
        fx -= bias;
        if ( fx > 0 ) {
            if ( dlabs[k] > 0 )
                ++ntp;
            else
                ++nfp;
        } else {
            if ( dlabs[k] > 0 )
                ++ntn;
            else
                ++nfn;
        }
    }
    in.close();
    analyze ( ntp, nfp, ntn, nfp );
        timer.stop();
        cerr << "prediction time  = " << timer.elapsedTime() << endl;
        return;
}


} // end namespace
#endif /* SVMPREDICTOR_H_ */
