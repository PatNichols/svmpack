/*
 * SVMKernelEvaluator.h
 *
 *  Created on: Jul 8, 2010
 *      Author: d3p708
 */

#ifndef SVMKERNELEVALUATOR_H_
#define SVMKERNELEVALUATOR_H_
#include "SVMOptions.h"
#include "svm_utils.h"
#include "svm_traits.h"
#include <cmath>
namespace svmpack
{

template<class svm_real>
class SVMKernelEvaluator
{
public:
    SVMKernelEvaluator ( const SVMOptions<svm_real>& options ) :
            vecs ( options.getVectorsPtr() ),
            scal ( 0x0 ),
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
        try {
            scal = new svm_real[nvecs];
        } catch ( exception& e ) {
            cerr << "allocation failed for scal factor in svm_kernel_matrix\n";
            cerr << "nvecs = " << nvecs << endl;
            exit ( EXIT_FAILURE );
        }
        const svm_real one ( 1 );
        if ( ktype == 2 ) {
            c1 = -c1;
        }
        if ( options.scaleKernel() ) {
            const svm_real * v1 = vecs;
            const svm_real tau = svmpack::svm_traits<svm_real>::tau();
            switch ( ktype ) {
            case 0:
                for ( size_t k = 0; k < nvecs; ++k, v1 += nfeat ) {
                    svm_real sk = eval0 ( v1, v1 );
                    if ( sk > tau ) {
                        scal[k] = one / sqrt ( sk );
                    } else {
                        scal[k] = one;
                    }
                }
                break;
            case 1:
                for ( size_t k = 0; k < nvecs; ++k, v1 += nfeat ) {
                    svm_real sk = eval1 ( v1, v1 );
                    if ( sk > tau ) {
                        scal[k] = one / sqrt ( sk );
                    } else {
                        scal[k] = one;
                    }
                }
                break;
            case 2:
                for ( size_t k = 0; k < nvecs; ++k ) {
                    scal[k] = one;
                }
                break;
            case 3:
                for ( size_t k = 0; k < nvecs; ++k, v1 += nfeat ) {
                    svm_real sk = eval3 ( v1, v1 );
                    if ( sk > tau ) {
                        scal[k] = one / sqrt ( sk );
                    } else {
                        scal[k] = one;
                    }
                }
                break;
            }
        } else {
            for ( size_t k = 0; k < nvecs; ++k ) {
                scal[k] = one;
            }
        }
    }
    ;

    ~SVMKernelEvaluator() {
        delete [] scal;
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

    inline svm_real eval ( const svm_real * __restrict__ v1, const size_t ix ) const throw () {
        const svm_real * __restrict__ v2 = vecs + ix * nfeat;
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

    inline void genRow ( svm_real * __restrict__ res, const size_t ix ) const throw () {
        const svm_real *const v1 = vecs + ix * nfeat;
        const svm_real s1 = scal[ix];
        const svm_real * v2 = vecs;
        size_t k;
        switch ( ktype ) {
        case 0:
#pragma omp parallel for private(k)         
            for ( k = 0; k < nvecs; ++k ) {
                res[k] = eval0 ( v1, v2 + k *nfeat ) * scal[k] * s1;
            }
            break;
        case 1:
#pragma omp parallel for private(k)         
            for ( size_t k = 0; k < nvecs; ++k ) {
                res[k] = eval1 ( v1, v2+k*nfeat ) * scal[k] * s1;
            }
            break;
        case 2:
#pragma omp parallel for private(k)         
            for ( size_t k = 0; k < nvecs; ++k ) {
                res[k] = eval2 ( v1, v2+k*nfeat );
            }
            break;
        case 3:
#pragma omp parallel for private(k)         
            for ( size_t k = 0; k < nvecs; ++k ) {
                res[k] = eval3 ( v1, v2+k*nfeat ) * scal[k] * s1;
            }
            break;
        }
        return;
    }
    ;
private:
    const svm_real *vecs;
    svm_real *scal;
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

}

#endif /* SVMKERNELEVALUATOR_H_ */
