/*
 * svm_utils.hpp
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef SVM_UTILS_HPP_
#define SVM_UTILS_HPP_
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "StringTokenizer.h"
#include <unistd.h>
using namespace std;

namespace svmpack
{

inline size_t getMemorySize()
{
    size_t npages = sysconf ( _SC_PHYS_PAGES );
    size_t pgsize = sysconf ( _SC_AVPHYS_PAGES );
    return npages * pgsize;
};

inline void analyze ( size_t ntp, size_t nfp, size_t ntn, size_t nfn )
{
    size_t nerr = nfp + ntn;
    size_t ncorr = nfn + ntp;
    size_t n = nerr + ncorr;
    size_t np = nfp + ntp;
    size_t nn = nfn + ntn;
    size_t nt = ntp + ntn;
    size_t nf = nfn + nfp;
    cerr << "# predictions          = " << n << endl;
    cerr << "# true                 = " << nt << endl;
    cerr << "# false                = " << nf << endl;
    cerr << "# positive             = " << np << endl;
    cerr << "# negative             = " << nn << endl;
    cerr << "# true-positive        = " << ntp << endl;
    cerr << "# true-negative        = " << ntn << endl;
    cerr << "# false-positive       = " << nfp << endl;
    cerr << "# false-negative       = " << nfn << endl;
    cerr << "accuracy               = " << ( double ( ncorr ) / double ( n ) ) << endl;
    cerr << "# errors               = " << nerr << endl;
    cerr << "precision              = " << ( double ( ntp ) / double ( np ) ) << endl;
    cerr << "recall                 = " << ( double ( ntp ) / double ( nt ) ) << endl;
    cerr << "sensitivity(true acc)  = " << ( double ( ntp ) / double ( nt ) ) << endl;
    cerr << "specificity(false acc) = " << ( double ( nfn ) / double ( nf ) ) << endl;
    cerr << "positive pred rate     = " << ( double ( ntp ) / double ( np ) ) << endl;
    cerr << "negative pred rate     = " << ( double ( nfn ) / double ( nn ) ) << endl;
    cerr << "False Positive rate    = " << ( double ( nfp ) / double ( nf ) ) << endl;
    cerr << "False Negative rate    = " << ( double ( ntn ) / double ( nt ) ) << endl;
    double sens = double ( ntp ) / double ( nt );
    double spec = double ( nfn ) / double ( nf );
    cerr << "Likelihood ratio pos   = " << ( sens / ( 1 - spec ) ) << endl;
    cerr << "Likelihood ratio neg   = " << ( ( 1 - sens ) / spec ) << endl;
    double prec = double ( ntp ) / double ( np );
    double recall = double ( ntp ) / double ( nt );
    double f = 2 * ( prec * recall ) / ( prec + recall );
    cerr << "F measure              = " << f << endl;
    double mcn = double ( ntp ) * double ( nfn ) - double ( nfp ) * double ( ntn );
    double mcd = sqrt ( double ( np ) * double ( nn ) * double ( nt ) * double ( nf ) );
    if ( mcd > 1.e-14 ) {
        mcn = mcn / mcd;
    } else {
        mcn = 1;
    }
    cerr << "Matthews's Correlation = " << mcn << endl;
}

template<class svm_real> inline svm_real powi ( svm_real x, size_t m )
{
    svm_real y;
    switch ( m ) {
    case 0:
        return 1;
    case 1:
        return x;
    case 2:
        return ( x * x );
    case 3:
        return ( x * x * x );
    case 4:
        x = x * x;
        return ( x * x );
    default:
        break;
    }
    y = ( m % 2 ) ? x : 1;
    while ( m >>= 1 ) {
        x *= x;
        if ( m % 2 )
            y *= x;
    }
    return y;
}

//////////////////////////////////////////////////////
// convert a libsvm data file to a tdo data file
//////////////////////////////////////////////////////
inline void libsvm2tdo ( const char *filename )
{
    string outfilename = string ( filename ) + string ( ".tdo" );
    ifstream in ( filename );
    ofstream out ( outfilename.c_str() );
    if ( !in ) {
        cerr << "could not open data file " << filename << "\n";
        exit ( EXIT_FAILURE );
    }
    if ( !out ) {
        cerr << "could not open output file " << outfilename << "\n";
        exit ( EXIT_FAILURE );
    }
    size_t nfeat ( 0 ), nvecs ( 0 );
    StringTokenizer toker ( "s", " :\n\t\0" );
    string sline;
    while ( getline ( in, sline ) ) {
        toker.resetString ( sline );
        if ( toker.hasMoreTokens() ) {
            string lab = toker.nextToken();
            ++nvecs;
            while ( toker.hasMoreTokens() ) {
                size_t indx = toker.nextElement<size_t> ();
                if ( toker.hasMoreTokens() ) {
                    string sval = toker.nextToken();
                } else {
                    cerr
                        << " no matching value pair for index in input file : "
                        << indx << " on line " << ( nvecs ) << endl;
                }
                if ( indx > nfeat )
                    nfeat = indx;
            }
        } else {
            break;
        }
    }
    in.clear();
    cerr << "read " << nvecs << " vector w/ " << nfeat << " features\n";
    in.seekg ( 0, ios::beg );
    double *vecs;
    double *yalf;
    try {
        vecs = new double[nvecs * nfeat];
        yalf = new double[nvecs];
    } catch ( exception& e ) {
        cerr << "error allocating memory in readDataFile\n";
        cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
        exit ( EXIT_FAILURE );
    }
    memset ( vecs, 0, sizeof ( double ) *nvecs * nfeat );
    size_t ivec ( 0 );
    while ( getline ( in, sline ) ) {
        toker.resetString ( sline );
        if ( toker.hasMoreTokens() ) {
            int il = toker.nextElement<int> ();
            if ( il > 0 ) {
                yalf[ivec] = double(1);
            } else {
                yalf[ivec] = double(-1);
            }
            while ( toker.hasMoreTokens() ) {
                size_t indx = toker.nextElement<size_t> ();
                --indx;
                if ( toker.hasMoreTokens() ) {
                    vecs[ivec * nfeat + indx] = toker.nextElement<double> ();
                } else {
                    cerr
                        << " no matching value pair for index in input file : "
                        << indx << " on line " << ( nvecs ) << endl;
                }
            }
            ++ivec;
        } else {
            break;
        }
    }
    in.close();
    int itmp;
    itmp = static_cast< size_t > ( nvecs );
    out.write ( ( char* ) &itmp, sizeof ( int ) );
    itmp = static_cast< size_t > ( nfeat );
    out.write ( ( char* ) &itmp, sizeof ( int ) );
    out.write ( ( char* ) yalf, sizeof ( double ) *nvecs );
    out.write ( ( char* ) vecs, sizeof ( double ) *nvecs * nfeat );
    out.close();
};

template <class real_t>
inline real_t fmax ( real_t x, real_t y ) throw()
{
    return (x>y)?x:y;
};
template <class real_t>
inline real_t fmin ( real_t x, real_t y ) throw()
{
    return (x<y)?x:y;
};

template < class T> inline T half() throw()
{
    return static_cast<T> ( 0.5 );
};
template <> inline float half<float>() throw()
{
    return 0.5f;
};
template <> inline double half<double>() throw()
{
    return 0.5;
};
template <> inline long double half<long double>() throw()
{
    return 0.5L;
};

}

#endif /* SVM_UTILS_HPP_ */
