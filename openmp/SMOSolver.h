/*
 * SMOSolver.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef SMOSOLVER_H_
#define SMOSOLVER_H_
#include "svm_utils.h"
#include "SVMOptions.h"
#include "SVMKernelMatrix.h"
#include "svm_traits.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <omp.h>
using namespace std;

namespace svmpack
{

template<class svm_real>
class SMOSolver
{
public:
    SMOSolver ( const SVMOptions<svm_real>& options_in ) :
            kmatrix ( options_in ),
            kmax(0x0),kmin(0x0),
            alpha ( 0x0 ), grad ( 0x0 ),
            dy ( options_in.getYAlphaPtr() ),
            y ( 0x0 ), status ( 0x0 ),
            nvecs ( options_in.getNVectors() ),
            maxits ( options_in.getMaxIterations() ),
            fsum ( 0 ), bias ( 0 ), gap ( 0 ),
            eps ( options_in.getEpsilon() ),
            cost ( options_in.getCost() ),
            options ( options_in ) {
        try {
            kmax=new svm_real*[1];
            kmin=new svm_real*[1];
            alpha = new svm_real[nvecs];
            grad = new svm_real[nvecs];
            y = new int[nvecs];
            status = new int[nvecs];
        } catch ( exception& e ) {
            cerr << "could not allocate arrays for smo solver\n";
            exit ( EXIT_FAILURE );
        }
        const svm_real zero ( 0 );
#ifdef _CHECK_DY_
        for ( size_t k = 0; k < nvecs; ++k ) {
            svm_real t=fabs(dy[k])-svm_real(1);
            if (t>eps) {
                cerr << " dy is out of whack ! " << dy[k] << " "<< k << endl;
            }
        }
#endif
        memset(alpha,0,sizeof(svm_real)*nvecs);
        memcpy(grad,dy,sizeof(svm_real)*nvecs);
        for ( size_t k=0; k<nvecs; ++k) {
            status[k] = -1;
            y[k] = ( dy[k] > zero ) ? 1:-1;
        }
    }
    ;

    ~SMOSolver() {
        delete[] status;
        delete[] y;
        delete[] grad;
        delete[] alpha;
        delete[] kmin;
        delete[] kmax;
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
        register size_t k;
        const svm_real zero ( 0 );
        register svm_real asum ( 0 ), csum ( 0 );
        register int nfree ( 0 );
        fsum = bias = svm_real ( 0 );
#pragma omp parellel for  private(k) reduction(+,asum) reduction(+,fsum)         
        for ( k = 0; k < nvecs; ++k ) {
            asum += alpha[k];
            fsum += alpha[k] * grad[k] * dy[k];
        }
        fsum = ( fsum + asum ) / svm_real ( 2 );
#pragma omp parellel for  private(k) reduction(+,bias) reduction(+,nfree)         
        for ( k = 0; k < nvecs; ++k ) {
            if ( status[k] != 0 )
                continue;
            bias += grad[k] ;
            ++nfree;
        }
        bias = -bias;
        if ( nfree )
            bias /= svm_real ( nfree );
#pragma omp parellel for  private(k) reduction(+,csum)         
        for ( k = 0; k < nvecs; ++k ) {
            csum += svmpack::fmax<svm_real> ( zero, ( dy[k] * ( grad[k] + bias ) ) );
        }
        csum *= cost;

        cerr << "csum = " <<csum<<endl;
        cerr << "asum = " <<asum<<endl;
        gap = ( csum + asum - fsum - fsum ) / ( svm_real ( 1 ) + asum + csum - fsum );
    };

    void takeStep ( const int imax, const int imin ) throw() {
        const svm_real TAU = svm_traits<svm_real>::tau();
        svm_real a1, a2, L, H, ai, aj;
        int st1, st2;
        register size_t k;
        const svm_real zero ( 0 );
        step_return = 0;
        int s = y[imax] * y[imin];
        svm_real ds = svm_real ( s );
        kmatrix.getRows ( kmax, kmin, imax, imin );
        const svm_real *qmax= *kmax;
        const svm_real *qmin= *kmin;
        a1 = ai = alpha[imax];
        a2 = aj = alpha[imin];
        const svm_real gam = a2 + ds * a1;
        if ( s > 0 ) {
            L = svmpack::fmax<svm_real> ( zero, ( gam - cost ) );
            H = svmpack::fmin<svm_real> ( cost, gam );
        } else {
            L = svmpack::fmax<svm_real> ( zero, ( gam ) );
            H = svmpack::fmin<svm_real> ( cost, ( gam + cost ) );
        }
        const svm_real qc = svmpack::fmax<svm_real> (
                                (qmax[imax] + qmin[imin] - qmax[imin] - qmax[imin]), TAU );
        a2 += ( ( grad[imin] - grad[imax] ) * dy[imin] ) / qc;
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
        if ( fabs ( a2 - aj ) > TAU ) {
            alpha[imax] = a1;
            alpha[imin] = a2;
            status[imax] = st1;
            status[imin] = st2;
            svm_real da1 = ( dy[imax] ) * ( a1 - ai );
            svm_real da2 = ( dy[imin] ) * ( a2 - aj );
#pragma omp parellel for  private(k) 
            for ( k = 0; k < nvecs; ++k ) {
                grad[k] -= ( da1 * qmax[k] + da2 * qmin[k] );
            }
            step_return = 1;
        }
    };


    void update() noexcept;

    ;
    void train() throw () {
        svm_real diff = 0;
        svm_real fold = 0;
        svm_stopwatch timer;
        timer.start();
        for ( size_t iter = 0; iter < maxits; ++iter ) {
            for ( size_t k = 0; k < nvecs; ++k ) {
                update();
                if ( step_return != 1 ) break;
            }
            findGap();
            diff = fsum - fold;
            fold = fsum;
            cerr << "iteration      = " << iter << endl;
            cerr << "obj. function  = " << fsum << endl;
            cerr << "diff. obj fun  = " << diff << endl;
            cerr << "gap            = " << gap << endl;
            cerr << "bias           = " << bias << endl << endl;
            if ( step_return != 1 ) {
                cerr << " converged! no more feasible step\n";
                break;
            }
            if ( gap < eps ) {
                cerr << " converged! gap is within tolerance\n";
                break;
            }
        }
        timer.stop();
        cerr << "training time             = " << timer.elapsedTime() << " seconds\n";
    };

    void outputModelFile() throw () {
        size_t nsv = 0;
        size_t nbnd = 0;
#pragma omp parellel for  private(k) reduction(+,nsv) reduction(+,nbnd)         
        for ( size_t k = 0; k < nvecs; ++k ) {
            if ( status[k] >= 0 ) {
                ++nsv;
                if ( status[k] > 0 ) ++nbnd;
            }
        }
        cerr << "# training vectors     = " << nvecs << endl;
        cerr << "# support vectors      = " << nsv << endl;
        cerr << "# bound support vectors= " << nbnd << endl;
        ofstream out ( options.getModelFileName().c_str() );
        if ( !out ) {
            cerr << "could not open file " << options.getModelFileName() << " to output model \n";
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
            if ( status[k] < 0 ) continue;
            dtmp = svm_real ( y[k] ) * kmatrix.getScaleFactor ( k ) * alpha[k];
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
                for ( size_t k = 0; k < ( nfeat ); ++k ) {
                    dtmp = static_cast<double> ( vp[k+j*nfeat] );
                    out.write ( ( char* ) &dtmp, sizeof ( double ) );
                }
            }
        }
        out.close();
        ofstream out2 ( options.getOutputFileName().c_str() );
        if ( !out2 ) {
            cerr << "could not open file " << options.getOutputFileName() << " to for output of scores \n";
            exit ( EXIT_FAILURE );
        }
        size_t ntp ( 0 ), nfp ( 0 ), ntn ( 0 ), nfn ( 0 );
        for ( size_t k = 0; k < nvecs; ++k ) {
            svm_real fx = ( dy[k] - grad[k] - bias );
            if ( fx > 0 ) {
                if ( y[k] > 0 ) ++ntp;
                else ++nfp;
            } else {
                if ( y[k] > 0 ) ++ntn;
                else ++nfn;
            }
            dtmp = static_cast<double> ( dy[k] );
            out2.write ( ( char* ) &dtmp, sizeof ( double ) );
            dtmp = static_cast<double> ( fx );
            out2.write ( ( char* ) &dtmp, sizeof ( double ) );
        }
        out2.close();
        analyze ( ntp, nfp, ntn, nfn );
    };
private:
    SVMKernelMatrix<svm_real> kmatrix;
    svm_real **kmax;
    svm_real **kmin;
    svm_real *alpha;
    svm_real *grad;
    const svm_real *dy;
    int *y;
    int *status;
    int step_return;
    size_t nvecs, maxits;
    svm_real fsum;
    svm_real bias;
    svm_real gap;
    svm_real eps, cost;
    const SVMOptions<svm_real>& options;
};


struct DoubleIndexPair {
    double value;
    std::int64_t index;
    
    DoubleIndexPair():value(0),index(-1) {}
    
    DoubleIndexPair(const double& d,const std::int64_t& i):value(d),index(i) {}
    
    DoubleIndexPair& Max(const DoubleIndexPair& x) noexcept {
        if (x.value > value) {
            value = x.value;
            index = x.index;
        }
        return *this;
    }
    DoubleIndexPair& Min(const DoubleIndexPair& x) noexcept {
        if (x.value < value) {
            value = x.value;
            index = x.index;
        }
        return *this;
    }    
    DoubleIndexPair& Max(const double& x, std::int64_t i) noexcept {
        if (x > value) {
            value = x;
            index = i;
        }
        return *this;
    }
    DoubleIndexPair& Min(const double& x,std::int64_t i) noexcept {
        if (x < value) {
            value = x;
            index = i;
        }
        return *this;
    }    
};

struct FloatIndexPair {
    float value;
    std::int64_t index;
    
    FloatIndexPair():value(0),index(-1) {}
    
    FloatIndexPair(const double& d,const std::int64_t& i):value(d),index(i) {}
    
    FloatIndexPair& Max(const FloatIndexPair& x) noexcept {
        if (x.value > value) {
            value = x.value;
            index = x.index;
        }
        return *this;
    }
    FloatIndexPair& Min(const FloatIndexPair& x) noexcept {
        if (x.value < value) {
            value = x.value;
            index = x.index;
        }
        return *this;
    }    
    FloatIndexPair& Max(const float& x, std::int64_t i) noexcept {
        if (x > value) {
            value = x;
            index = i;
        }
        return *this;
    }
    FloatIndexPair& Min(const float& x,std::int64_t i) noexcept {
        if (x < value) {
            value = x;
            index = i;
        }
        return *this;
    }    
};


template < typename T >
inline void SMOSolver<T>::update() noexcept {
        step_return = 0;
        register T gmin = svm_traits<T>::huge();
        register T gmax = -gmin;
        register int imax = -1;
        register int imin = -1;
        int max_list[128];
        int min_list[128];
        T gmax_l[128];
        T gmin_l[128];
        T m_gmax,m_gmin;
        T m_imax,m_imin; 
        std::size_t k;   
#pragma omp parallel private(m_imax,m_imin,m_gmax,m_gmin) shared(gmax_l,gmin_l,max_list,min_list)
   {
        int tid = omp_get_thread_num();
        int nth = omp_get_num_threads();
        int m_imax = imax;
        int m_imin = imin;
        T m_gmax = gmax;
        T m_gmin = gmin;         
#pragma omp parallel for shared(nvecs) private(k)  
        for ( k = 0; k < nvecs; ++k ) {
            register int ys = y[k] * status[k];
            register T gk = grad[k];
            if ( ys != 1 && gk > m_gmax ) {
                m_gmax = gk;
                m_imax = k;
            }
            if ( ys != -1 && gk < m_gmin ) {
                m_gmin = gk;
                m_imin = k;
            }
        }
        max_list[tid]=m_imax;
        min_list[tid]=m_imin;
        gmax_l[tid]=m_gmax;
        gmin_l[tid]=m_gmin;
#pragma omp single 
        {
            gmax = gmax_l[0];
            imax = max_list[0];
            for (std::size_t k=1;k<nth;++k) {
                T gmx = gmax_l[k];
                if (gmx > gmax) {
                    gmax=gmx;
                    imax=max_list[k];
                }
            }
        }
#pragma omp single
        {
            gmin = gmin_l[0];
            imin = min_list[0];
            for (std::size_t k=1;k<nth;++k) {
                T gmx = gmin_l[k];
                if (gmx < gmin) {
                    gmin=gmx;
                    imin=min_list[k];
                }
            }

        }
   }
        if ( imax != -1 && imin != -1 && imax != imin && ( ( gmax - gmin ) > eps ) ) {
            SMOSolver<T>::takeStep ( imax, imin );
        }
    }
    ;


template <>
inline void SMOSolver<double>::update() noexcept {
        step_return = 0;
        const std::int64_t lneg1 = std::int64_t{1}; 

        DoubleIndexPair mymax(-2.e300,-1LL);
        DoubleIndexPair mymin( 2.e300,-1LL);
#pragma omp declare reduction( MyDMax: DoubleIndexPair: omp_out.Max(omp_in)) \
 initializer ( omp_priv=DoubleIndexPair(-2.e300,-1LL) )

#pragma omp declare reduction( MyDMin: DoubleIndexPair: omp_out.Min(omp_in)) \
 initializer ( omp_priv=DoubleIndexPair( 2.e300,-1LL) )

       std::size_t k;
#pragma omp parallel for shared(nvecs) private(k) reduction(MyDMax:mymax) reduction(MyDMin:mymin)
        for ( k = 0; k < nvecs; ++k ) {
            auto ys = y[k] * status[k];
            if ( ys != 1 ) mymax.Max(grad[k],k);
            if ( ys != -1) mymin.Min(grad[k],k);
        }
        if ( mymax.index != -1 && mymin.index != -1 && mymax.index != mymin.index && ( ( mymax.value - mymin.value ) > eps ) ) {
            SMOSolver<double>::takeStep ( mymax.index, mymin.index );
        }
    }
    ;


template <>
inline void SMOSolver<float>::update() noexcept {
        step_return = 0;
        const std::int64_t lneg1 = std::int64_t{1}; 

        FloatIndexPair mymax(-2.e30f,-1LL);
        FloatIndexPair mymin( 2.e30f,-1LL);
#pragma omp declare reduction( MyDMax: FloatIndexPair: omp_out.Max(omp_in)) \
 initializer ( omp_priv=FloatIndexPair(-2.e30f,-1LL) )

#pragma omp declare reduction( MyDMin: FloatIndexPair: omp_out.Min(omp_in)) \
 initializer ( omp_priv=FloatIndexPair( 2.e30f,-1LL) )

       std::int64_t k;
#pragma omp parallel for private(k) shared(nvecs) reduction(MyDMax:mymax) reduction(MyDMin:mymin)
        for ( k = 0; k < nvecs; ++k ) {
            auto ys = y[k] * status[k];
            if ( ys != 1 ) mymax.Max(grad[k],k);
            if ( ys != -1) mymin.Min(grad[k],k);
        }
        if ( mymax.index != -1 && mymin.index != -1 && mymax.index != mymin.index && ( ( mymax.value - mymin.value ) > eps ) ) {
            SMOSolver<float>::takeStep ( mymax.index, mymin.index );
        }
    }
    ;


}

#endif /* SMOSOLVER_H_ */