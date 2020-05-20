/*
 * SVMKernelMatrix.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef SVMKERNELMATRIX_H_
#define SVMKERNELMATRIX_H_
#include "SVMOptions.h"
#include "SVMKernelEvaluator.h"
#include "svm_stopwatch.h"
#include <cstdlib>
using namespace std;

namespace svmpack
{

template<class svm_real>
class SVMKernelMatrixThread
{
public:
    SVMKernelMatrixThread ( const SVMOptions<svm_real>& options,
                            svm_real *kmat_in, size_t thread_id ) :
            keval ( new SVMKernelEvaluator<svm_real> ( options ) ),
            kmat ( kmat_in ),
            nvecs ( options.getNVectors() ),
            nth ( options.getNThreads() ),
            tid ( thread_id ),
            off1 ( 0 ), off2 ( 0 ) {
        size_t bsz = nvecs / nth;
        size_t xsz = nvecs % nth;
        if ( tid < xsz ) {
            off1 = tid * ( bsz + 1 );
            off2 = off1 + bsz + 1;
        } else {
            off1 = bsz * tid + xsz;
            off2 = off1 + bsz;
        }
    };

    ~SVMKernelMatrixThread() {
        delete keval;
    }

    void run() throw() {
#ifdef _USE_BUFFER_
        size_t BUFFER_SIZE = (128L * 1048576L);
        size_t bsize = BUFFER_SIZE/nvecs/sizeof(svm_real);
        while (bsize<1) {
            BUFFER_SIZE <<=1;
            bsize = BUFFER_SIZE/nvecs/sizeof(svm_real);
        }
        if (bsize > (off2-off1)) bsize=koff2-koff1; 
        const size_t asize = bsize * nvecs;
        svm_real *buffer= new svm_real[asize];
        svm_real *brow=buffer;
        register size_t cnt=0;
        for (register size_t k= off1; k<off2; ++k) {
            keval->genRow( brow, k);
            brow+=nvecs;
            ++cnt;
            if (cnt==bsize) {
                memcpy( kmat + (k-cnt+1)*nvecs, buffer, sizeof(svm_real)*cnt*nvecs);
                brow=buffer;
                cnt=0;
            }
        }
        if (cnt) {
            memcpyh(kmat+(koff2-cnt)*nvecs,buffer,sizeof(svm_real)*cnt*nvecs);
        }
        delete [] buffer;
#else
        svm_real *kmat_k = kmat + off1 * nvecs;
        for ( size_t k = off1; k < off2; ++k, kmat_k += nvecs ) keval->genRow ( kmat_k, k );
#endif
    };

private:
    const SVMKernelEvaluator<svm_real> *keval;
    svm_real *kmat;
    size_t nvecs, nth, tid;
    size_t off1, off2;
};

namespace kthx
{
template < class svm_real > static void * krun ( void *arg )
{
    SVMKernelMatrixThread<svm_real> *t =
        reinterpret_cast< SVMKernelMatrixThread< svm_real >* > ( arg );
    t->run();
    return 0x0;
}
}

template<class svm_real>
class SVMKernelMatrix
{
public:
    SVMKernelMatrix ( const SVMOptions<svm_real>& options ) :
            kfun ( options ),
            data ( 0x0 ), indx ( 0x0 ),
            last ( 0 ), nvecs ( options.getNVectors() ),
            cache_size ( options.getCacheSize() ) {
        svm_stopwatch timer;
        timer.start();
        if ( cache_size < 2 ) cache_size = 0;
        try {
            size_t asize = ( cache_size == 0 ) ? nvecs : cache_size;
            asize *= nvecs;
            data = new svm_real[asize];
        } catch ( exception& e ) {
            cerr << "allocation failed for scal factor in svm_kernel_matrix\n";
            cerr << "nvecs = " << nvecs << endl;
            cerr << "cache_size" << cache_size << endl;
            exit ( EXIT_FAILURE );
        }
        if ( cache_size != 0 ) {
            try {
                indx = new int[cache_size];
            } catch ( exception& e ) {
                cerr << "allocation failed in " << __FILE__ << " ON LINE "
                     << __LINE__ << endl;
                exit ( EXIT_FAILURE );
            }
            for ( size_t k = 0; k < cache_size; ++k ) {
                indx[k] = -1;
            }
            indx[0] = 0;
            last = 0;
            kfun.genRow ( data, 0 );
            cerr << " cache initialized!\n";
        } else {
            size_t nth = options.getNThreads();
            if ( nth == 0 ) {
                std::size_t k;
                cerr << "serial kernel matrix\n";
                svm_real * datak = data;
#pragma omp paralell for private(k)
                for ( k = 0; k < nvecs; ++k) {
                    kfun.genRow ( ( data + k * nvecs ), k );
                }
            } else {
                cerr << "threaded kernel matrix\n";
                pthread_t *pth = new pthread_t[nth];
                SVMKernelMatrixThread<svm_real> **t =
                    new SVMKernelMatrixThread<svm_real>*[nth];
                for ( size_t k = 0; k < nth; ++k )
                    t[k] = new SVMKernelMatrixThread<svm_real> ( options, data, k );
                for ( size_t k = 0; k < nth; ++k )
                    pthread_create ( ( pth + k ), 0x0, kthx::krun<svm_real>, ( void* ) t[k] );
                for ( size_t k = 0; k < nth; ++k )
                    pthread_join ( pth[k], 0x0 );
                for ( size_t k = 0; k < nth; ++k )
                    delete t[k];
                delete [] t;
                delete [] pth;
            }
            cerr << " kernel matrix computed!\n";
        }
        timer.stop();
        cerr << "initialized kernel matrix in " << timer.elapsedTime() << " seconds\n";
    }
    ;

    ~SVMKernelMatrix() {
        if ( indx ) delete [] indx;
        delete[] data;
    }

    inline svm_real getScaleFactor ( size_t ix ) const throw () {
        return kfun.getScaleFactor ( ix );
    }
    ;

    inline void getRows ( svm_real **  qmax,
                          svm_real **  qmin,
                          int imax, int imin ) throw () {
        if ( cache_size == 0 ) {
            *qmax = data + imax * nvecs;
            *qmin = data + imin * nvecs;
            return;
        }
        size_t max_found = 0;
        size_t min_found = 0;
        size_t k;
        for ( k = 0; k < cache_size; ++k ) {
            int ik = indx[k];
            if ( ik == imax ) max_found = k + 1;
            if ( ik == imin ) min_found = k + 1;
        }
        if ( max_found ) {
            --max_found;
            *qmax = data + ( max_found ) * nvecs;
        } else {
            last = ( last + 1 ) % cache_size;
            svm_real *dp = data + last * nvecs;
            kfun.genRow ( dp, imax );
            *qmax = dp;
            indx[last] = imax;
        }
        if ( min_found ) {
            --min_found;
            *qmin = data + ( min_found ) * nvecs;
        } else {
            last = ( last + 1 ) % cache_size;
            svm_real *dp = data + last * nvecs;
            kfun.genRow ( dp, imin );
            *qmin = dp;
            indx[last] = imin;
        }
    }
    ;

    const svm_real *getScaleFactorPtr() const throw() {
        return kfun.getScaleFactorPtr();
    };

    const svm_real *getKernelMatrixPtr() const throw() {
        if ( cache_size == 0 ) return data;
        cerr << "can not return kernel matrix ptr in SVMKernelMatrix::getKernelMatrixPtr()";
        cerr << "\n we are using a cache!\n";
        exit ( EXIT_FAILURE );
    };
private:
    SVMKernelEvaluator< svm_real> kfun;
    svm_real *data;
    int *indx;
    size_t last, nvecs, cache_size;
};

}

#endif /* SVMKERNELMATRIX_H_ */
