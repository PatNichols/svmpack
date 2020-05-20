/*
 * svm_traits.hpp
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef SVM_TRAITS_HPP_
#define SVM_TRAITS_HPP_
#include <cstdlib>
#include <limits>
using namespace std;

namespace svmpack
{

template <class T> struct svm_traits {
    static inline  T eps() throw() {
        return numeric_limits<T>::epsilon() * 100;
    };
    static inline  T tau() throw() {
        return numeric_limits<T>::epsilon * 10;
    }
    static inline  T huge() throw() {
        return numeric_limits<T>::huge();
    };
};

template <> struct svm_traits<double> {
    static inline  double eps() throw() {
        return 2.e-12;
    };
    static inline  double tau() throw() {
        return 2.e-14;
    }
    static inline  double huge() throw() {
        return 2.e100;
    };
};

template <> struct svm_traits<float> {
    static inline  float eps() throw() {
        return 2.e-6f;
    };
    static inline  float tau() throw() {
        return 2.e-7f;
    }
    static inline  float huge() throw() {
        return 2.e30f;
    };
};


}

#endif /* SVM_TRAITS_HPP_ */
