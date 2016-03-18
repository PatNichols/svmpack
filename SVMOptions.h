/*
 * SVMOptions.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef SVMOPTIONS_H_
#define SVMOPTIONS_H_
#include "svm_traits.h"
#include "svm_utils.h"
#include "ProgramOptions.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;

namespace svmpack
{

template<class svm_real>
class SVMOptions
{
public:
    SVMOptions ( int argc, char **argv );

    ~SVMOptions() {
        delete[] vecs;
        delete[] yalf;
    }
    ;

    void readDataFile() throw ();
    void readLIBSVMFile() throw ();
    void readModelFile() throw ();

    svm_real getKernelCof1() const throw () {
        return cof1;
    }
    ;

    svm_real getKernelCof2() const throw () {
        return cof2;
    }
    ;

    string getDataFileName() const throw () {
        return datafile;
    }
    ;

    svm_real getEpsilon() const throw () {
        return eps;
    }
    ;

    size_t getKernelType() const throw () {
        return ktype;
    }
    ;

    size_t getKernelPower() const throw () {
        return kpow;
    }
    ;

    size_t getMaxIterations() const throw () {
        return maxits;
    }
    ;

    string getModelFileName() const throw () {
        return modelfile;
    }
    ;

    size_t getNFeatures() const throw () {
        return nfeat;
    }
    ;

    size_t getNVectors() const throw () {
        return nvecs;
    }
    ;

    size_t getCacheSize() const throw () {
        return cache_size;
    }
    ;

    string getOutputFileName() const throw () {
        return outfile;
    }
    ;

    bool scaleKernel() const throw () {
        return scale;
    }
    ;

    const svm_real *getVectorsPtr() const throw () {
        return vecs;
    }
    ;

    const svm_real *getYAlphaPtr() const throw () {
        return yalf;
    }
    ;

    size_t getNThreads() const throw () {
        return nproc;
    }
    ;

    size_t getTask() const throw () {
        return task;
    }
    ;

    svm_real getCost() const throw () {
        return cost;
    }
    ;

    svm_real getBias() const throw () {
        return bias;
    }
    ;

    ostream& writeToStream ( ostream& os ) const throw() {
        size_t ntrue = 0;
        for ( size_t k = 0; k < nvecs; ++k ) if ( yalf[k] > 0 ) ++ntrue;
        size_t nfalse = nvecs - ntrue;
        os << "task                 = ";
        switch ( task ) {
        case 0:
            cerr << "train" << endl;
            break;
        case 1:
            cerr << "predict" << endl;
            break;
        case 2:
            cerr << "classify" << endl;
            break;
        }
        os << "# vectors            = " << nvecs << endl;
        os << "# features           = " << nfeat << endl;
        os << "# true               = " << ntrue << endl;
        os << "# false              = " << nfalse << endl;
        os << "# threads            = " << nproc << endl;
        os << "data file            = " << datafile << endl;
        os << "model file           = " << modelfile << endl;
        os << "output file          = " << outfile << endl;
        os << "kernel type          = " << ktype << endl;
        os << "kernel cof1          = " << cof1 << endl;
        os << "kernel cof2          = " << cof2 << endl;
        os << "kernel power         = " << kpow << endl;
        os << "kernel scale         = " << scale << endl;
        if ( task == 0 ) {
            os << "kernel cache size    = " << cache_size << endl;
            os << "solver eps           = " << eps << endl;
            os << "solver max iterations= " << maxits << endl;
            os << "solver cost          = " << cost << endl;
        }
        if ( task == 1 ) {
            os << "bias                 = " << bias << endl;
        }
        os << "sizeof(SVM_REAL)     = " << sizeof ( svm_real ) << endl;
        return os;
    };
private:
    size_t nfeat, nvecs, cache_size;
    size_t maxits, nproc, task;
    size_t ktype, kpow;
    svm_real * vecs;
    svm_real * yalf;
    svm_real cof1, cof2;
    svm_real eps, cost, bias;
    bool scale;
    string modelfile;
    string datafile;
    string outfile;

};

template <class svm_real> ostream& operator<< ( ostream& os,
        const SVMOptions<svm_real>& opts )
{
    return opts.writeToStream ( os );
}

template<class svm_real> inline SVMOptions<svm_real>::SVMOptions ( int argc, char **argv )
{
    ProgramOptions opts;
    ostringstream ostr;
    ostr << svmpack::svm_traits<svm_real>::eps();
    string epsstr = ostr.str();
    string
    kdes =
        "kernel function type\n\
\t\t\t\t0) dot product\n\
\t\t\t\t1) polynomial k(x,y)=pow( c1*dot(x,y)+c2,pow)\n\
\t\t\t\t2) Gaussian k(x,y)= exp( -c1 * dot (x-y,x-y) )\n\
\t\t\t\t4) Sigmoid  k(x,y)= tanh (c1*dot(x,y)+c2)";
    opts.addOption ( "kernel_type", kdes.c_str(), "2" );
    opts.addOption ( "kernel_power", " pow for kernel type 2", "2" );
    opts.addOption ( "kernel_scale",
                     " scale kernel elements so  kernel matrix diagonal is unity",
                     "true" );
    opts.addOption ( "kernel_cof1", " c1 in the above formulas for the kernel function" );
    opts.addOption (
        "kernel_cof2",
        " c2 in the above formulas for the kernel function (kernel type=1)",
        "1" );
    opts.addOption ( "eps", "training convergence parameter", epsstr.c_str() );
    opts.addOption ( "cost",
                     "soft margin penalty for training (infinity is hard margin case)",
                     "1" );
    opts.addOption (
        "cache_size",
        "number of kernel matrix rows to store in cache(<2 precompute the kernel)",
        "0" );
    opts.addOption ( "max_iterations", "max number of training iterations" );
    opts.addOption ( "nthreads", "number of threads to use(0 is serial run)", "0" );
    opts.addOption ( "task", "task to perform (train,predict,classify)", "train" );
    opts.addOption ( "data", "name of datafile", "svm.in" );
    opts.addOption ( "model", "name of model file", "svm.model" );
    opts.addOption ( "out", "name of file to output labels and scores", "svm.out" );
    opts.addOption ( "config", "name of file to read for config options" );
    opts.parseCommandLine ( argc, argv );
    if ( opts.hasValue ( "config" ) ) {
        string cfile = opts.getValue<string> ( "config" );
        opts.parseConfigFile ( cfile.c_str() );
    }
    datafile = opts.getValue<string> ( "data" );
    modelfile = opts.getValue<string> ( "model" );
    outfile = opts.getValue<string> ( "out" );
    eps = opts.getValue<svm_real> ( "eps" );
    string task_str = opts.getValue<string> ( "task" );
    if ( task_str[0] == 't' ) {
        task = 0;
        if ( datafile.find ( ".tdo" ) != string::npos ) {
            cerr << "reading TDO file = " << datafile << endl;
            readDataFile();
        } else {
            cerr << "reading libsvm file = " << datafile << endl;
            readLIBSVMFile();
        }
        if ( nfeat == 0 || nvecs == 0 ) {
            cerr << "no input data was read!\n";
            cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << "\n";
            exit ( EXIT_FAILURE );
        }
        ktype = opts.getValue<size_t> ( "kernel_type" );
        kpow = opts.getValue<size_t> ( "kernel_power" );
        scale = opts.getValue<bool> ( "kernel_scale" );
        if ( opts.hasValue ( "kernel_cof1" ) ) {
            cof1 = opts.getValue<svm_real> ( "kernel_cof1" );
        } else {
            if ( ktype == 2 )
                cof1 = svm_real ( 10 ) / svm_real ( nfeat );
            else
                cof1 = 1;
        }
        cof2 = opts.getValue<svm_real> ( "kernel_cof2" );
        cost = opts.getValue<svm_real> ( "cost" );
        if ( opts.hasValue ( "max_iterations" ) ) {
            maxits = opts.getValue<size_t> ( "max_iterations" );
        } else {
            maxits = nvecs;
        }
        nproc = opts.getValue<size_t> ( "nthreads" );
        cache_size = opts.getValue<size_t> ( "cache_size" );
        size_t msize = svmpack::getMemorySize();
        if ( cache_size == 0 ) {
            size_t ksize = sizeof ( svm_real ) * nvecs * nvecs;
            if ( ksize > msize ) {
                cerr
                    << " KERNEL MATRIX TOO LARGE TO BE STORES IN AVAILABLE MEMORY!!!!\n";
                cerr << " using cache !\n";
                while ( ( ksize >>= 1 ) > msize )
                    ;
                cache_size = ksize / sizeof ( svm_real ) / nvecs;
                if ( cache_size < 2 ) {
                    cerr << " data too large for cache even!!!\n";
                    cerr << " nvecs = " << nvecs << "\n";
                    exit ( EXIT_FAILURE );
                }
            }
        }
        if ( cache_size == 1 )
            cache_size = 2;
        if ( cache_size > 1 ) {
            if ( cache_size > nvecs )
                cache_size = nvecs;
            size_t csize = cache_size * sizeof ( svm_real ) * nvecs;
            if ( csize > msize ) {
                cerr << "Given cache size " << cache_size
                     << " requires too much memory!\n";
                while ( ( csize >>= 1 ) > msize )
                    ;
                cache_size = csize / sizeof ( svm_real ) * nvecs;
                if ( cache_size < 2 ) {
                    cerr << "cannot hold 2 kernel rows in memory\n";
                    cerr << "cache_size= " << cache_size << endl;
                    exit ( EXIT_FAILURE );
                }
            }
        }
    }
    if ( task_str[0] == 'p' ) {
        task = 1;
        nproc = opts.getValue<size_t> ( "nthreads" );        
        readModelFile();
    }
    if ( task_str[0] == 'c' ) {
        task = 2;
        nproc = opts.getValue<size_t> ( "nthreads" );
        readModelFile();
    }
}

template<class svm_real> inline void svmpack::SVMOptions<svm_real>::readDataFile() throw()
{
    ifstream in ( datafile.c_str() );
    if ( !in ) {
        cerr << "could not open data file " << datafile << "\n";
        exit ( EXIT_FAILURE );
    }
    int itmp;
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    nvecs = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    nfeat = static_cast<size_t> ( itmp );
    try {
        size_t vsize = nvecs * nfeat;
        vecs = new svm_real[vsize];
        yalf = new svm_real[nvecs];
    } catch ( ... ) {
        cerr << "error allocating memory in readDataFile\n";
        cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
        exit ( EXIT_FAILURE );
    }
    if ( sizeof ( svm_real ) == sizeof ( double ) ) {
        in.read ( ( char* ) yalf, sizeof ( double ) * nvecs );
        in.read ( ( char* ) vecs, sizeof ( double ) * nvecs * nfeat );
    } else {
        double *dv = new double[nvecs];
        if ( !dv ) {
            cerr << "error allocating memory in readDataFile\n";
            cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
            exit ( EXIT_FAILURE );
        }
        in.read ( ( char* ) dv, sizeof ( double ) * nvecs );
        for ( size_t k = 0; k < nvecs; ++k )
            yalf[k] = static_cast<svm_real> ( dv[k] );
        svm_real *vptr = vecs;
        for ( size_t k = 0; k < nfeat; ++k ) {
            in.read ( ( char* ) dv, sizeof ( double ) *nvecs );
            for ( size_t l = 0; l < nvecs; ++l, ++vptr )
                ( *vptr ) = static_cast<svm_real> ( dv[l] );
        }
        delete[] dv;
    }
    in.close();
}

template<class svm_real> inline void svmpack::SVMOptions<svm_real>::readLIBSVMFile() throw()
{
    ifstream in ( datafile.c_str() );
    if ( !in ) {
        cerr << "could not open data file " << datafile << "\n";
        exit ( EXIT_FAILURE );
    }
    nfeat = nvecs = 0;
    StringTokenizer toker ( "s", " :\n\t\0" );
    string sline;
    while ( getline ( in, sline ) ) {
        toker.resetString ( sline );
        if ( toker.hasMoreTokens() ) {
            string lab = toker.nextToken();
            ++nvecs;
            while ( toker.hasMoreTokens() ) {
                size_t indx = toker.nextElement<size_t>();
                if ( toker.hasMoreTokens() ) {
                    string sval = toker.nextToken();
                } else {
                    cerr << " no matching value pair for index in input file : " << datafile <<
                         " on line " << ( nvecs ) << endl;
                }
                if ( indx > nfeat ) nfeat = indx;
            }
        } else {
            break;
        }
    }
    in.clear();
    cerr << "read " << nvecs << " vector w/ " << nfeat << " features\n";
    in.seekg ( 0, ios::beg );
    try {
        size_t vsize = nvecs;
        vsize *= nfeat;
        vecs = new svm_real[vsize];
        yalf = new svm_real[nvecs];
    } catch ( exception& e ) {
        cerr << "error allocating memory in readDataFile\n";
        cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
        exit ( EXIT_FAILURE );
    }
    memset ( vecs, 0, sizeof ( svm_real ) *nvecs * nfeat );
    for ( size_t k = 0; k < nvecs; ++k ) {
        getline ( in, sline );
        toker.resetString ( sline );
        if ( toker.hasMoreTokens() ) {
            int il = toker.nextElement<int>();
            if ( il > 0 ) {
                yalf[k] = svm_real(1);
            } else {
                yalf[k] = svm_real(-1);
            }
            while ( toker.hasMoreTokens() ) {
                size_t indx = toker.nextElement<size_t>();
                --indx;
                if ( toker.hasMoreTokens() ) {
                    vecs[ ( k*nfeat+indx ) ] = toker.nextElement<svm_real>();
                } else {
                    cerr << " no matching value pair for index in input file : " << datafile <<
                         " on line " << ( nvecs ) << endl;
                }
            }
        } else {
            break;
        }
    }
    in.close();
}

template<class svm_real> inline void svmpack::SVMOptions<svm_real>::readModelFile() throw()
{
    ifstream in ( modelfile.c_str() );
    if ( !in ) {
        cerr << "could not open model file " << modelfile << "\n";
        exit ( EXIT_FAILURE );
    }
    int itmp;
    double dtmp;
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    nvecs = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    nfeat = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    ktype = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    kpow = static_cast<size_t> ( itmp );
    in.read ( ( char* ) &dtmp, sizeof ( double ) );
    cof1 = static_cast<svm_real> ( dtmp );
    in.read ( ( char* ) &dtmp, sizeof ( double ) );
    cof2 = static_cast<svm_real> ( dtmp );
    in.read ( ( char* ) &dtmp, sizeof ( double ) );
    bias = static_cast<svm_real> ( dtmp );
    in.read ( ( char* ) &itmp, sizeof ( int ) );
    if ( itmp == 0 )
        scale = false;
    else
        scale = true;
    try {
        size_t vsize = nvecs * nfeat;
        vecs = new svm_real[vsize];
        yalf = new svm_real[nvecs];
    } catch ( exception& e ) {
        cerr << "error allocating memory in readModelFile\n";
        cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
        exit ( EXIT_FAILURE );
    }
    if ( sizeof ( svm_real ) == sizeof ( double ) ) {
        in.read ( ( char* ) yalf, sizeof ( double ) * nvecs );
        in.read ( ( char* ) vecs, sizeof ( double ) * nvecs * nfeat );
    } else {
        double *dv = new double[nvecs];
        if ( !dv ) {
            cerr << "error allocating memory in readModelFile\n";
            cerr << "nvecs = " << nvecs << " nfeat = " << nfeat << endl;
            exit ( EXIT_FAILURE );
        }
        in.read ( ( char* ) dv, sizeof ( double ) * nvecs );
        for ( size_t k = 0 ; k < nvecs; ++k ) {
            yalf[k] = static_cast< svm_real > ( dv[k] );
        }
        for ( size_t k = 0; k < nfeat; ++k ) {
            in.read ( ( char* ) dv, sizeof ( double ) * nvecs );
            for ( size_t l = 0; l < nvecs; ++l )
                vecs[ ( k * nvecs + l ) ] = static_cast<svm_real> ( dv[l] );
        }
        delete[] dv;
    }
    in.close();
}

}
#endif /* SVMOPTIONS_H_ */
