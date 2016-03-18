/*
 * svm_stopwatch.h
 *
 *  Created on: Jul 8, 2010
 *      Author: d3p708
 */

#ifndef SVM_STOPWATCH_H_
#define SVM_STOPWATCH_H_
#define _USE_RT_
#ifdef _USE_RT_
#include <ctime>
namespace svmpack
{
struct svm_stopwatch {
    double acc;
    struct timespec ts;
    struct timespec tf;

    svm_stopwatch() : acc ( 0 ) {};

    ~svm_stopwatch() {};

    void start() throw() {
        clock_gettime ( CLOCK_MONOTONIC, &ts );
    };

    void stop() throw() {
        clock_gettime ( CLOCK_MONOTONIC, &tf );
        acc += diffTime ( tf, ts );
    };

    void clear() throw() {
        acc = 0;
    };

    double elapsedTime() const throw() {
        return acc;
    };

    double diffTime ( struct timespec& tf, struct timespec& ts ) {
        double diff = tf.tv_sec - ts.tv_sec;
        diff += 1.e-9 * ( tf.tv_nsec - ts.tv_nsec );
        return diff;
    };

};
}
#else
#include <sys/time.h>
namespace svmpack
{
struct svm_stopwatch {
    double acc;
    struct timeval ts, tf;

    svm_stopwatch() : acc ( 0 ) {};

    ~svm_stopwatch() {};

    void start() throw() {
        gettimeofday ( &ts, 0x0 );
    };

    void stop() throw() {
        gettimeofday ( &tf, 0x0 );
        acc += diffTime ( tf, ts );
    };

    void clear() throw() {
        acc = 0;
    };

    double elapsedTime() const throw() {
        return acc;
    };

    static inline double diffTime ( struct timeval& tf, struct timeval& ts ) {
        double diff = tf.tv_sec - ts.tv_sec;
        diff += 1.e-6 * ( tf.tv_usec - ts.tv_usec );
        return diff;
    };
};
}
#endif

#endif /* SVM_STOPWATCH_H_ */
