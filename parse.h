/*
 * parse.hpp
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef PARSE_HPP_
#define PARSE_HPP_
#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>
#include <cfloat>
#include <cerrno>
using namespace std;
namespace svmpack
{

template<class T> inline T parse ( const char *str ) throw ()
{
    istringstream in ( string ( str ) );
    T x;
    in >> x;
    if ( in )
        return x;
    cerr << "error parse " << str << "\n";
    exit ( EXIT_FAILURE );
}

template<> inline float parse<float> ( const char *str ) throw ()
{
    char *end_ptr;
    errno = 0;
    float val = strtof ( str, &end_ptr );
    if ( end_ptr == str ) {
        cerr << "no digits found for float\n";
        exit ( EXIT_FAILURE );
    }
    if ( errno ) {
        cerr << "error in reading float\n";
        exit ( EXIT_FAILURE );
    }
    return val;
}

template<> inline double parse<double> ( const char *str ) throw ()
{
    char *end_ptr;
    errno = 0;
    double val = strtod ( str, &end_ptr );
    if ( end_ptr == str ) {
        cerr << "no digits found for double\n";
        exit ( EXIT_FAILURE );
    }
    if ( errno ) {
        cerr << "error in reading double\n";
        exit ( EXIT_FAILURE );
    }
    return val;
}

template<> inline long parse<long> ( const char *str ) throw ()
{
    char *end_ptr;
    errno = 0;
    long val = strtol ( str, &end_ptr, 10 );
    if ( end_ptr == str ) {
        cerr << "no digits found for long\n";
        exit ( EXIT_FAILURE );
    }
    if ( ( errno == ERANGE && ( val == LONG_MAX || val == LONG_MIN ) ) || ( errno != 0
            && val == 0 ) ) {
        cerr << "error in reading long\n";
        exit ( EXIT_FAILURE );
    }
    return val;
}

template<> inline unsigned long parse<unsigned long> ( const char *str ) throw ()
{
    char *end_ptr;
    errno = 0;
    unsigned long val = strtoul ( str, &end_ptr, 10 );
    if ( end_ptr == str ) {
        cerr << "no digits found for ulong\n";
        exit ( EXIT_FAILURE );
    }
    if ( ( errno == ERANGE && ( val == ULONG_MAX ) ) || ( errno != 0 && val == 0 ) ) {
        cerr << "error in reading ulong\n";
        exit ( EXIT_FAILURE );
    }
    return val;
}

template<> inline int parse<int> ( const char *str ) throw ()
{
    char *end_ptr;
    errno = 0;
    long val = strtol ( str, &end_ptr, 10 );
    if ( end_ptr == str ) {
        cerr << "no digits found for int\n";
        exit ( EXIT_FAILURE );
    }
    if ( ( errno == ERANGE && ( val == LONG_MAX || val == LONG_MIN ) ) || ( errno != 0
            && val == 0 ) ) {
        cerr << "error in reading int\n";
        exit ( EXIT_FAILURE );
    }
    if ( val > INT_MAX || val < INT_MIN ) {
        cerr << "cannot convert to int\n";
        exit ( EXIT_FAILURE );
    }
    return static_cast<int> ( val );
}

template<> inline unsigned int parse<unsigned int> ( const char *str ) throw ()
{
    char *end_ptr;
    errno = 0;
    unsigned long val = strtoul ( str, &end_ptr, 10 );
    if ( end_ptr == str ) {
        cerr << "no digits found for uint\n";
        exit ( EXIT_FAILURE );
    }
    if ( ( errno == ERANGE && ( val == ULONG_MAX ) ) || ( errno != 0 && val == 0 ) ) {
        cerr << "error in reading uint\n";
        exit ( EXIT_FAILURE );
    }
    if ( val > ULONG_MAX ) {
        cerr << "cannot convert to uint\n";
        exit ( EXIT_FAILURE );
    }
    return static_cast<unsigned int> ( val );
}

template <> inline string parse<string> ( const char *str ) throw()
{
    return string ( str );
}

template <> inline bool parse<bool> ( const char *str ) throw()
{
    char ch = str[0];
    if ( ch == 'T' || ch == 't' || ch == '1' ) return true;
    if ( ch == 'F' || ch == 'f' || ch == '0' ) return false;
    if ( string ( str ).find_first_not_of ( "0" ) == string::npos ) {
        return false;
    }
    return true;
}

}

#endif /* PARSE_HPP_ */
