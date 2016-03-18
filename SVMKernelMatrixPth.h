/*
 * SVMKernelMatrixPth.h
 *
 *  Created on: Jul 9, 2010
 *      Author: d3p708
 */
#ifndef SVMKERNELMATRIXPTH_H_
#define SVMKERNELMATRIXPTH_H_
#include <pthread.h>
#include "SVMBarrierPth.h"
#include <cmath>
#include <cstdlib>
#include "SVMOptions.h"
#include "svm_traits.h"
#include "svm_stopwatch.h"

using namespace std;
namespace svmpack
{

template<class svm_real>
class SVMKernelEvaluatorPth
{
public:
    SVMKernelEvaluatorPth ( const SVMOptions<svm_real>& options, size_t thread_id ) :
            vecs ( options.getVectorsPtr() ),
            scal ( 0x0 ),
            c1 ( options.getKernelCof1() ),
            c2 ( options.getKernelCof2() ),
            nvecs ( options.getNVectors() ),
            nfeat ( options.getNFeatures() ),
            ktype ( options.getKernelType() ),
            kpow ( options.getKernelPower() ),
            nth ( options.getNThreads() ),
            tid ( thread_id ),
            off1 ( 0 ), off2 ( 0 ) {
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
        size_t bsz = nvecs / nth;
        size_t xsz = nvecs % nth;
        if ( tid < xsz ) {
            off1 = bsz * tid + tid;
            off2 = off1 + bsz + 1;
        } else {
            off1 = bsz * tid + xsz;
            off2 = off1 + bsz;
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

    ~SVMKernelEvaluatorPth() {
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
        const svm_real *__restrict__ const v1 = vecs + ix * nfeat;
        const svm_real s1 = scal[ix];
        const svm_real *__restrict__  v2 = vecs + off1 * nfeat;
        switch ( ktype ) {
        case 0:
            for ( size_t k = off1; k < off2; ++k, v2 += nfeat ) {
                res[k] = eval0 ( v1, v2 ) * scal[k] * s1;
            }
            break;
        case 1:
            for ( size_t k = off1; k < off2; ++k, v2 += nfeat ) {
                res[k] = eval1 ( v1, v2 ) * scal[k] * s1;
            }
            break;
        case 2:
            for ( size_t k = off1; k < off2; ++k, v2 += nfeat ) {
                res[k] = eval2 ( v1, v2 );
            }
            break;
        case 3:
            for ( size_t k = off1; k < off2; ++k, v2 += nfeat ) {
                res[k] = eval3 ( v1, v2 ) * scal[k] * s1;
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
    size_t nth, tid, off1, off2;

    inline svm_real eval0 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
        svm_real s ( 0 ) __attribute__((aligned(32)));
        for ( size_t k = 0; k < nfeat; ++k )
            s += v1[k] * v2[k];
        return s;
    }
    ;

    inline svm_real eval1 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
        svm_real s ( 0 );
        for (  size_t k = 0; k < nfeat; ++k )
            s += v1[k] * v2[k];
        s=c1*s+c2;
        if (kpow==2) {
            return s*s;
        }
        return svmpack::powi<svm_real> (s, kpow );
    }
    ;

    inline svm_real eval3 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
         svm_real s ( 0 );
        for (  size_t k = 0; k < nfeat; ++k )
            s += v1[k] * v2[k];
        return tanh ( ( c1 * s + c2 ) );
    }
    ;

    inline svm_real eval2 ( const svm_real * __restrict__ v1,
                            const svm_real * __restrict__ v2 ) const throw () {
         svm_real t, s ( 0 );
        for (  size_t k = 0; k < nfeat; ++k ) {
            t = v1[k] - v2[k];
            s += t * t;
        }
        s *= c1;
        return exp ( s );
    }
    ;

};



template<class svm_real>
class SVMKernelMatrixPth
{
public:

    SVMKernelMatrixPth ( const SVMOptions<svm_real>& options ) :
            kfun ( options, 0 ),
            barrier ( static_cast<int>(options.getNThreads()),
                      mem_fun ( &SVMKernelMatrixPth<svm_real>::updateCache ) ),
            data ( 0x0 ), indx ( 0x0 ),
            is_found ( 0x0 ), is_valid ( 0x0 ),
            max_curr ( 0x0 ), min_curr ( 0x0 ),
            coff1 ( 0x0 ), coff2 ( 0x0 ), last ( 0x0 ),
            cache_size ( options.getCacheSize() ),
            nvecs ( options.getNVectors() ) {
        if ( cache_size < 2 ) {
            cerr << "cache size less than 2 in SVMKernelMatrixPth!\n";
            exit ( EXIT_FAILURE );
        }
        svm_stopwatch timer;
        timer.start();
        try {
            data = new svm_real[ ( nvecs * cache_size ) ];
            indx = new int[cache_size];
            is_found = new size_t[2];
            is_valid = new size_t[2];
            max_curr = new size_t ( 0 );
            min_curr = new size_t ( 0 );
            last = new size_t ( 0 );
        } catch ( ... ) {
            cerr << "allocation failed for scal factor in svm_kernel_matrix\n";
            cerr << "nvecs = " << nvecs << endl;
            exit ( EXIT_FAILURE );
        }
        size_t nth = options.getNThreads();
        size_t bsz = cache_size / nth;
        size_t xsz = cache_size % nth;
        if ( xsz ) {
            coff1 =  0;
            coff2 =  bsz + 1;
        } else {
            coff1 =  0;
            coff2 = bsz;
        }
        for ( size_t k = 0; k < cache_size; ++k ) {
            indx[k] = -1;
        }
        *last = 0;
        kfun.genRow ( data, 0 );
        is_found[0] = is_found[1] = 0;
        timer.stop();
        cerr << "initialized kernel matrix in " << timer.elapsedTime()
             << " seconds\n";
    }
    ;

    SVMKernelMatrixPth ( const SVMOptions<svm_real>& options,
                         SVMKernelMatrixPth<svm_real>& c, size_t thread_id ) :
            kfun ( options, thread_id ),
            barrier ( c.barrier, static_cast<int>(thread_id) ),
            data ( c.data ),
            indx ( c.indx ),
            is_found ( c.is_found ),
            is_valid ( c.is_valid ),
            max_curr ( c.max_curr ),
            min_curr ( c.min_curr ),
            coff1 ( 0x0 ),
            coff2 ( 0x0 ),
            last ( c.last ),
            cache_size ( c.cache_size ),
            nvecs ( c.nvecs ) {
        svm_stopwatch timer;
        timer.start();
        size_t nth = options.getNThreads();
        size_t bsz = cache_size / nth;
        size_t xsz = cache_size % nth;
        if ( thread_id < xsz ) {
            coff1 = bsz * thread_id + thread_id;
            coff2 = coff1 + bsz + 1;
        } else {
            coff1 = bsz * thread_id + xsz;
            coff2 = coff1 + bsz;
        }
        kfun.genRow ( data, 0 );
        timer.stop();
        cerr << "initialized kernel matrix in " << timer.elapsedTime()
             << " seconds\n";
    }
    ;

    ~SVMKernelMatrixPth() {
        if ( coff1 == 0 ) {
            delete last;
            delete min_curr;
            delete max_curr;
            delete[] is_found;
            delete[] indx;
            delete[] data;
        }
    };

    inline svm_real getScaleFactor ( size_t ix ) const throw () {
        return kfun.getScaleFactor ( ix );
    };


    inline const svm_real *getScaleFactorPtr() const throw () {
        return kfun.getScaleFactorPtr();
    }
    ;

    void updateCache() throw() {
        int key=is_found[1]-1;
        if ( is_found[0] ) {
            *max_curr = is_found[0] - 1;
            is_valid[0] = 1;
        } else {
            do {
                *last = ( *last + 1 ) % cache_size;
            } while (key == *last);
            *max_curr = *last;
            indx[*last] = imax_tot;
            is_valid[0] = 0;
        }
        key= *max_curr;
        if ( is_found[1] ) {
            *min_curr = is_found[1] - 1;
            is_valid[1] = 1;
        } else {
            do {
                *last = ( *last + 1 ) % cache_size;
            } while (*last==key);
            *min_curr = *last;
            indx[*last] = imin_tot;
            is_valid[1] = 0;
        }
        is_found[0] = 0;
        is_found[1] = 0;
    };

    inline void getRows ( svm_real ** qmax,
                          svm_real ** qmin,
                          int imax, int imin ) throw () {
        imax_tot = imax;
        imin_tot = imin;
        if ( imax != -1 && imin != -1 && imax != imin ) {
            for (  size_t k = coff1; k < coff2; ++k ) {
                 int ik = indx[k];
                if ( ik == imax ) is_found[0] = k + 1;
                if ( ik == imin ) is_found[1] = k + 1;
            }
            barrier.reduce ( this );
            svm_real *dmax = data + ( *max_curr ) * nvecs;
            svm_real *dmin = data + ( *min_curr ) * nvecs;
            if ( !is_valid[0] ) kfun.genRow ( dmax, imax );
            if ( !is_valid[1] ) kfun.genRow ( dmin, imin );
            *qmax = dmax;
            *qmin = dmin;
        }
    }
    ;
private:
    SVMKernelEvaluatorPth<svm_real> kfun;
    svmpack::CyclicBarrierPth< SVMKernelMatrixPth<svm_real> > barrier;
    svm_real *data;
    int  *indx;
    size_t *is_found;
    size_t *is_valid;
    size_t *max_curr;
    size_t *min_curr;
    size_t imax_tot;
    size_t imin_tot;
    size_t coff1, coff2, *last;
    size_t cache_size, nvecs;
};

}

#endif /* SVMKERNELMATRIXPTH_H_ */
