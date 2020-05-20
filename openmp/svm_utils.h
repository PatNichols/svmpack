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
#include <unistd.h>
#include <vector>
#include "parse.h"
using namespace std;

namespace svmpack
{

inline size_t getMemorySize()
{
    return size_t(16ULL*1024ULL*1048576ULL);
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
    cerr << "# accuracy             = " << ( double ( ncorr ) / double ( n ) ) << endl;
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
    size_t nfeat=0;
    size_t nvecs=0;
    string sline;
    vector<string> tokens;
    while ( getline ( in, sline ) ) {
        explodeString(sline," :\n\0",tokens);
        size_t s = tokens.size();
        ++nvecs;
        if (s%2==0) {
            std::cerr << "parse error on line " << nvecs << "\n";
            std::cerr << "no label or bad index:value pair\n";
            exit(EXIT_FAILURE);
        }
        size_t last_index = stoull(tokens[s-2]);
        if (last_index > nfeat) nfeat = last_index;
    }
    in.clear();
    cerr << "read " << nvecs << " vector w/ " << nfeat << " features\n";
    in.seekg(0);
    double *vecs;
    double *yalf;
    size_t vsize = nvecs;
    vsize *=nfeat;
    try {
        vecs = new double[vsize];
        yalf = new double[nvecs];
    } catch ( exception& e ) {
        cerr << "error allocating memory in readDataFile\n";
        cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
        exit ( EXIT_FAILURE );
    }
    for (size_t i=0;i<vsize;++i) {
        vecs[i]=double(0);
    }
    for ( size_t k = 0; k < nvecs; ++k ) {
        getline ( in, sline );
        explodeString(sline," :\n",tokens);
        std::vector<std::string>::iterator iter = tokens.begin();
        std::vector<std::string>::iterator iend = tokens.end(); 
        int lab = stoi(*iter);
        if (lab>0) yalf[k]=1;
        else yalf[k]=-1;
        ++iter;
        while (iter!=iend) {
            size_t indx = stoull(*iter);
            --indx;
            ++iter;
            double d = stod(*iter);
            ++iter;
            vecs[k*nfeat+indx] = static_cast<double>(d);
        }
    }
    in.close();
    int itmp;
    itmp = static_cast< size_t > ( nvecs );
    out.write ( ( char* ) &itmp, sizeof ( int ) );
    itmp = static_cast< size_t > ( nfeat );
    out.write ( ( char* ) &itmp, sizeof ( int ) );
    out.write ( ( char* ) yalf, sizeof ( double ) *nvecs );
    out.write ( ( char* ) vecs, sizeof ( double ) *vsize );
    out.close();
};

//////////////////////////////////////////////////////
// convert a tdo data file to a libsvm data file
//////////////////////////////////////////////////////
inline void tdo2libsvm ( const char *filename )
{
    string infile(filename);
    size_t pos = infile.find(".");
    if (pos==string::npos) {
        std::cerr << "tdo file is not properly named! should be name.tdo \n";
        std::cerr << "name is " << infile << "\n";
        exit(EXIT_FAILURE);
    }
    string outfilename = infile.substr(0,pos);  
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
    size_t nvecs,nfeat;
    int itmp;
    itmp = static_cast< size_t > ( nvecs );

    in.read ( ( char* ) &itmp, sizeof ( int ) );
    nvecs = itmp;
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    nfeat = itmp;
    double *yalf = new double[nvecs];
    size_t vsize = nfeat * nvecs;
    double *vecs = new double[vsize];
    in.read ( ( char* ) yalf, sizeof ( double ) *nvecs );
    in.read ( ( char* ) vecs, sizeof ( double ) *nvecs * nfeat );
    in.close();

    for (size_t k=0;k<nvecs;++k) {
        if (yalf[k] > 0.0) {
            out << " 1 ";
        }else{
            out << " -1 "; 
        }
        double *vk = vecs + k * nfeat;
        for (size_t m=0;m<nfeat;++m) {
            if (fabs(vk[m]) > 1.e-14) {
                out << (m+1) << ":" << vk[m] << " ";
            }
        }
        out <<"\n";
    }
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
