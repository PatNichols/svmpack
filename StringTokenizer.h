/*
 * StringTokenizer.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef STRINGTOKENIZER_H_
#define STRINGTOKENIZER_H_
#include <cstring>
#include <string>
#include "parse.h"
namespace svmpack
{

class StringTokenizer
{
#define _DEFAULT_DELIMITER_ " \t\r\n\v"
public:
    StringTokenizer ( string str_in, char *delimiter_in = 0x0 ) :
            str ( str_in ), delimiter ( delimiter_in ), first ( 0 ), last ( 0 ) {
        if ( delimiter == 0x0 ) delimiter = _DEFAULT_DELIMITER_;
        tokenize();
    };

    ~StringTokenizer() {};

    bool hasMoreTokens() const throw() {
        return ( first != string::npos );
    };

    size_t countTokens() throw() {
        size_t f = first;
        size_t l = last;
        size_t ntoken = 0;
        while ( first != string::npos ) {
            ++ntoken;
            tokenize();
        }
        first = f;
        last = l;
        return ntoken;
    };

    template <class T> inline T nextElement() throw() {
        return parse<T> ( this->nextToken().c_str() );
    };

    string nextToken() throw() {
        size_t f = first;
        size_t l = last;
        tokenize();
        return str.substr ( f, ( l - f ) );
    };

    void resetString ( string& str_in ) throw() {
        str = str_in;
        first = 0;
        last = 0;
        tokenize();
    };

    void resetDelimiter ( char *delimiter_in ) throw() {
        delimiter = delimiter_in;
        if ( delimiter == 0 ) {
            delimiter = _DEFAULT_DELIMITER_ ;
        }
    };
private:
    string str;
    char *delimiter;
    size_t first, last;

    void tokenize() {
        first = str.find_first_not_of ( delimiter, last );
        last = str.find_first_of ( delimiter, first );
    };

};


}

#endif /* STRINGTOKENIZER_H_ */
