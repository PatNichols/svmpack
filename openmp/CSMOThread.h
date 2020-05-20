/*
 * CSMOThread.h
 *
 *  Created on: Jul 13, 2010
 *      Author: d3p708
 */

#ifndef CSMOTHREAD_H_
#define CSMOTHREAD_H_
#include <pthread.h>
#include <cstdlib>
#include "SVMBarrierPth.h"
#include "SVMOptions.h"
#include <iostream>
#include <iomanip>

using namespace std;

namespace svmpack
{

template<class svm_real>
struct CSMOThread {

    CSMOThread ( const SVMOptions<svm_real>& options,
                 const svm_real& cost_in,
                 const svm_real& eps_in,
                 svm_real *alpha_in,
                 svm_real *grad_in,
                 const svm_real *dy_in,
                 const int *y_in,
                 int *status_in,
                 size_t nVectors,
                 size_t maxIterations,
                 int nthreads ) :
            kmatrix ( options ),
            cost ( cost_in ), eps ( eps_in ), alpha ( alpha_in ), grad ( grad_in ), dy ( dy_in ),
            y ( y_in ), status ( status_in ), nvecs ( nVectors ),
            maxits ( maxIterations ), nth ( nthreads ), tid ( 0 ), off1 ( 0 ), off2 ( 0 ),
            step_return ( new size_t ( 0 ) ), old_alpha ( new svm_real[2] ),
            gmax_blk ( new svm_real[nthreads] ),
            gmin_blk ( new svm_real[nthreads] ),
            imax_blk ( new int[nthreads] ),
            imin_blk ( new int[nthreads] ),
            imax_tot ( new int ( -1 ) ),
            imin_tot ( new int ( -1 ) ),
            nsum_blk ( new int[nthreads] ),
            asum_blk ( new svm_real[nthreads] ),
            fsum_blk ( new svm_real[nthreads] ),
            bsum_blk ( new svm_real[nthreads] ),
            csum_blk ( new svm_real[nthreads] ),
            gap ( new svm_real ( 0 ) ),
            asum_tot ( new svm_real ( 0 ) ),
            fsum_tot ( new svm_real ( 0 ) ),
            bsum_tot ( new svm_real ( 0 ) ),
            kmax ( new svm_real*[1] ),
            kmin ( new svm_real*[1] ),
            maxmin ( nthreads, mem_fun ( &CSMOThread::reduceMaxMin ) ),
            step ( nthreads, mem_fun ( &CSMOThread::takeStep ) ),
            gap1 ( nthreads, mem_fun ( &CSMOThread::reduceBias ) ),
            gap2 ( nthreads, mem_fun ( &CSMOThread::reduceGap ) ) {
        size_t bsz = nvecs / nth;
        size_t xsz = nvecs % nth;
        if ( tid < xsz ) {
            off1 = bsz * tid + tid;
            off2 = off1 + bsz + 1;
        } else {
            off1 = bsz * tid + xsz;
            off2 = off1 + bsz;
        }
    }
    ;

    CSMOThread ( const SVMOptions<svm_real>& options, CSMOThread& t, int thread_id ) :
            kmatrix ( options, t.kmatrix, thread_id ),
            cost ( t.cost ), eps ( t.eps ), alpha ( t.alpha ), grad ( t.grad ), dy ( t.dy ),
            y ( t.y ), status ( t.status ), nvecs ( t.nvecs ),
            maxits ( t.maxits ), nth ( t.nth ), tid ( thread_id ), off1 ( 0 ), off2 ( 0 ),
            step_return ( t.step_return ), old_alpha ( t.old_alpha ),
            gmax_blk ( t.gmax_blk ), gmin_blk ( t.gmin_blk ),
            imax_blk ( t.imax_blk ), imin_blk ( t.imin_blk ),
            imax_tot ( t.imax_tot ), imin_tot ( t.imin_tot ),
            nsum_blk ( t.nsum_blk ), asum_blk ( t.asum_blk ),
            fsum_blk ( t.fsum_blk ), bsum_blk ( t.bsum_blk ),
            csum_blk ( t.csum_blk ), gap ( t.gap ),
            asum_tot ( t.asum_tot ),
            fsum_tot ( t.fsum_tot ), bsum_tot ( t.bsum_tot ),
            kmax ( t.kmax ),
            kmin ( t.kmin ),
            maxmin ( t.maxmin, thread_id ),
            step ( t.step, thread_id ),
            gap1 ( t.gap1, thread_id ),
            gap2 ( t.gap2, thread_id ) {
        size_t bsz = nvecs / nth;
        size_t xsz = nvecs % nth;
        if ( tid < xsz ) {
            off1 = bsz * tid + tid;
            off2 = off1 + bsz + 1;
        } else {
            off1 = bsz * tid + xsz;
            off2 = off1 + bsz;
        }
    }
    ;


    void takeStep() throw() {
        register svm_real a1, a2, ai, aj, L, H;
        register int st1, st2;
        const svm_real TAU = svm_traits<svm_real>::tau();
        const svm_real zero ( 0 );
        const int i = *imax_tot;
        const int j = *imin_tot;
        const svm_real * const qi = *kmax;
        const svm_real * const qj = *kmin;
        *step_return = 0;
        if ( i != j && i != -1 && j != -1 ) {
            int s = y[i] * y[j];
            svm_real ds ( s );
            a1 = ai = alpha[i];
            a2 = aj = alpha[j];
            old_alpha[0] = ai;
            old_alpha[1] = aj;
            const svm_real gam = a2 + ds * a1;
            if ( s > 0 ) {
                L = svmpack::fmax<svm_real> ( zero, ( gam - cost ) );
                H = svmpack::fmin<svm_real> ( cost, gam );
            } else {
                L = svmpack::fmax<svm_real> ( zero, ( gam ) );
                H = svmpack::fmin<svm_real> ( cost, ( gam + cost ) );
            }
            const svm_real qc = svmpack::fmax<svm_real> ( qi[i] + qj[j] - qi[j] - qi[j], TAU );
            a2 += ( ( grad[j] - grad[i] ) * dy[j] ) / qc;
            a2 = svmpack::fmin<svm_real> ( svmpack::fmax<svm_real> ( a2, L ), H );
            a1 += ds * ( aj - a2 );
            if ( a1 > TAU ) {
                if ( a1 < ( cost - TAU ) ) {
                    st1 = 0x0;
                } else {
                    a1 = cost;
                    a2 = gam - ds * cost;
                    st1 = 0x1;
                }
            } else {
                a1 = zero;
                a2 = gam;
                st1 = -0x1;
            }
            if ( a2 > TAU ) {
                if ( a2 < ( cost - TAU ) ) {
                    st2 = 0x0;
                } else {
                    st2 = 0x1;
                }
            } else {
                st2 = -0x1;
            }
            if ( fabs ( a2 - aj ) > eps ) {
                alpha[i] = a1;
                alpha[j] = a2;
                status[i] = st1;
                status[j] = st2;
                *step_return = 1;
            }
        }
    }
    ;

    void update() throw () {
        svm_real gmin = svm_traits<svm_real>::huge();
        svm_real gmax = -gmin;
        int imax = -1;
        int imin = -1;
        for ( size_t k = off1; k < off2; ++k ) {
            int ys = y[k] * status[k];
            if ( ys != 1 && grad[k] > gmax ) {
                gmax = grad[k];
                imax = k;
            }
        }
        for ( size_t k = off1; k < off2; ++k ) {
            int ys = y[k] * status[k];
            if ( ys != -1 && grad[k] < gmin ) {
                gmin = grad[k];
                imin = k;
            }
        }
        imax_blk[tid] = imax;
        imin_blk[tid] = imin;
        gmax_blk[tid] = gmax;
        gmin_blk[tid] = gmin;
        maxmin.reduce ( this );
        kmatrix.getRows ( kmax, kmin, *imax_tot, *imin_tot );
        step.reduce ( this );
        if ( *step_return == 1 ) {
            imax = *imax_tot;
            imin = *imin_tot;
            const svm_real * const qmax = *kmax;
            const svm_real * const qmin = *kmin;
            const svm_real da1 = ( dy[imax] ) * ( alpha[imax] - old_alpha[0] );
            const svm_real da2 = ( dy[imin] ) * ( alpha[imin] - old_alpha[1] );
            for ( register size_t k = off1; k < off2; ++k ) {
                grad[k] -= ( da1 * qmax[k] + da2 * qmin[k] );
            }
        }
    }
    ;


    void reduceMaxMin() throw() {
        register svm_real gmin = svm_traits<svm_real>::huge();
        register svm_real gmax = -gmin;
        register int imax = -1;
        register int imin = -1;
        register size_t k;
        for ( k = 0; k < nth; ++k ) {
            if ( gmax < gmax_blk[k] ) {
                gmax = gmax_blk[k];
                imax = imax_blk[k];
            }
            if ( gmin > gmin_blk[k] ) {
                gmin = gmin_blk[k];
                imin = imin_blk[k];
            }
        }
        *imax_tot = imax;
        *imin_tot = imin;
    };

    void reduceBias() throw() {
        register svm_real asum ( 0 ), bsum ( 0 );
        register int nsum ( 0 );
        for ( register size_t k = 0; k < nth; ++k ) {
            asum += asum_blk[k];
            bsum += bsum_blk[k];
            nsum += nsum_blk[k];
        }
        asum_tot[0] = asum;
        if ( nsum ) {
            bsum_tot[0] = bsum / nsum;
        } else {
            bsum_tot[0] = 0;
        }
    };

    void reduceGap() throw() {
        register svm_real csum ( 0 ), fsum ( 0 );
        for ( register size_t k = 0; k < nth; ++k ) {
            fsum += fsum_blk[k];
            csum += csum_blk[k];
        }
        svm_real tmp = asum_tot[0] + csum - fsum;
        gap[0] = ( tmp - fsum ) / ( tmp + 1 );
        fsum_tot[0] = fsum;
    };

    void findGap() throw () {
        register size_t k;
        register svm_real asum ( 0 ), csum ( 0 ), bsum ( 0 ), fsum ( 0 );
        register int nsum ( 0 );
        fsum = 0;
        bsum = 0;
        for ( k = off1; k < off2; ++k ) {
            asum += alpha[k];
            fsum += alpha[k] * grad[k] * dy[k];
        }
        fsum = ( fsum + asum ) * svmpack::half<svm_real>();
        for ( k = off1; k < off2; ++k ) {
            if ( status[k] )
                continue;
            bsum += grad[k];
            nsum += 1;
        }
        asum_blk[tid] = asum;
        bsum_blk[tid] = -bsum;
        fsum_blk[tid] = fsum;
        nsum_blk[tid] = nsum;
        gap1.reduce ( this );
        bsum = *bsum_tot;
        svm_real zero ( 0 );
        for ( k = off1; k < off2; ++k ) {
            csum += svmpack::fmax<svm_real> ( zero, ( ( grad[k] + bsum ) * dy[k] ) );
        }
        csum *= cost;
        csum_blk[tid] = csum;
        gap2.reduce ( this );
    }
    ;

    void run() throw() {
        svm_real diff ( 0 );
        svm_real fold ( 0 );
        for ( size_t iter = 0; iter < maxits; ++iter ) {
            for ( size_t k = 0; k < nvecs; ++k ) {
                update();
                if ( *step_return != 1 )
                    break;
            }
            findGap();
            diff = *fsum_tot - fold;
            fold = *fsum_tot;
            if ( tid == 0 ) {
                cerr << "iteration      = " << iter << endl;
                cerr << "obj. function  = " << *fsum_tot << endl;
                cerr << "diff. obj fun  = " << diff << endl;
                cerr << "gap            = " << *gap << endl;
                cerr << "bias           = " << *bsum_tot << endl << endl;
            }
            if ( *step_return != 1 ) {
                if ( !tid ) cerr << " converged! no more feasible step\n";
                break;
            }
            if ( *gap < eps ) {
                if ( !tid ) cerr << " converged! gap is within tolerance\n";
                break;
            }
        }
    }
    ;

    const svm_real* getScaleFactorPtr() const throw() {
        return kmatrix.getScaleFactorPtr();
    };

private:
    SVMKernelMatrixPth<svm_real> kmatrix;
    svm_real cost, eps;
    svm_real *alpha;
    svm_real *grad;
    const svm_real *dy;
    const int *y;
    int *status;
    size_t nvecs, maxits;
    int nth,tid;
    size_t off1,off2;
    size_t *step_return;
    svm_real *old_alpha;
    svm_real *gmax_blk;
    svm_real *gmin_blk;
    int *imax_blk;
    int *imin_blk;
    int *imax_tot;
    int *imin_tot;
    int *nsum_blk;
    svm_real *asum_blk;
    svm_real *fsum_blk;
    svm_real *bsum_blk;
    svm_real *csum_blk;
    svm_real *gap;
    svm_real *asum_tot;
    svm_real *fsum_tot;
    svm_real *bsum_tot;
    svm_real **kmax;
    svm_real **kmin;
    CyclicBarrierPth<CSMOThread<svm_real> > maxmin;
    CyclicBarrierPth<CSMOThread<svm_real> > step;
    CyclicBarrierPth<CSMOThread<svm_real> > gap1;
    CyclicBarrierPth<CSMOThread<svm_real> > gap2;
};

}

#endif /* CSMOTHREAD_H_ */
