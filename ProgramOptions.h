/*
 * ProgramOptions.h
 *
 *  Created on: Jul 7, 2010
 *      Author: d3p708
 */

#ifndef PROGRAMOPTIONS_H_
#define PROGRAMOPTIONS_H_
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "parse.h"
namespace svmpack
{

class ProgramOptions
{
    struct Option {
        Option ( const char *keyword, const char *description ) :
                key ( keyword ), des ( description ), val ( string ( "" ) ), stat ( 0 ) {
        }
        ;

        Option ( const char *keyword, const char *description, const char *value ) :
                key ( keyword ), des ( description ), val ( value ), stat ( -1 ) {
        }
        ;

        bool hasValue() const throw() {
            return ( stat != 0 );
        };

        bool matches ( const char *keyword ) const throw() {
            return ( key.compare ( keyword ) == 0 );
        }
        ;

        string getValue() const throw () {
            if ( stat ) {
                return string ( val );
            }
            cerr << " value requested but none was set in ProgramOption::Option::getValue()\n";
            exit ( EXIT_FAILURE );
        }
        ;

        void setValue ( const char *valueString ) throw () {
            if ( stat != 1 ) {
                val = valueString;
            }
            stat = 1;
        }
        ;

        ostream& writeToStream ( ostream& os ) const throw() {
            os << "-" << key << " " << des << "\n";
            if ( stat == -1 ) {
                os << "\t\t\t\t (default = " << val << ")\n";
            }
            return os;
        };

        string key;
        string des;
        string val;
        int stat;
    };

    vector<Option> options;
public:
    ProgramOptions() {
    }
    ;

    ~ProgramOptions() {
    }
    ;

    void parseConfigFile ( const char *filename ) throw() {
        ifstream in ( filename );
        if ( !in ) {
            cerr << "could not open the config file " << filename << "\n";
            exit ( EXIT_FAILURE );
        }
        while ( in ) {
            string s1, s2;
            in >> s1;
            if ( s1.size() == 0 )
                break;
            if ( !in.fail() ) {
                in >> s2;
                cerr << s1 << " " << s2 << endl;
                if ( !in.fail() ) setValue ( s1.c_str(), s2.c_str() );
            } else {
                cerr << "could not read the config file " << filename << "\n";
                exit ( EXIT_FAILURE );
            }
        }
        in.close();
    }
    ;

    void parseCommandLine ( int argc, char **argv ) throw() {
        for ( int k = 1; k < argc; k += 2 ) {
            if ( argv[k][0] == '-' ) {
                if ( strcmp ( argv[k], "-help" ) != 0 ) {
                    char *str = argv[k] + 1;
                    setValue ( str, argv[k + 1] );
                } else {
                    printHelp();
                }
            } else {
                cerr << "bad parse of command line options\n";
                cerr << "option " << argv[k] << endl;
                exit ( EXIT_FAILURE );
            }
        }
    }
    ;

    void addOption ( const char *keyword, const char *description ) {
        options.push_back ( Option ( keyword, description ) );
    }
    ;

    void addOption ( const char *keyword, const char *description,
                     const char *value ) {
        options.push_back ( Option ( keyword, description, value ) );
    }

    void setValue ( const char *keyword, const char *value ) throw() {
        vector<Option>::iterator iter = options.begin();
        vector<Option>::const_iterator iend = options.end();
        while ( iter != iend ) {
            if ( iter->matches ( keyword ) ) {
                iter->setValue ( value );
                return;
            }
            ++iter;
        }
        cerr << "could not find the option " << keyword << endl;
        exit ( EXIT_FAILURE );
    }
    ;

    template<class T> inline T getValue ( const char *keyword ) const throw() {
        vector<Option>::const_iterator iter = options.begin();
        vector<Option>::const_iterator iend = options.end();
        while ( iter != iend ) {
            if ( iter->matches ( keyword ) ) {
                string s = iter->getValue();
                return parse<T> ( s.c_str() );
            }
            ++iter;
        }
        cerr << "could not find the option " << keyword << endl;
        exit ( EXIT_FAILURE );
    }
    ;


    inline bool hasValue ( const char *keyword ) const throw() {
        vector<Option>::const_iterator iter = options.begin();
        vector<Option>::const_iterator iend = options.end();
        while ( iter != iend ) {
            if ( iter->matches ( keyword ) ) {
                return iter->hasValue();
            }
            ++iter;
        }
        cerr << "could not find the option " << keyword << endl;
        exit ( EXIT_FAILURE );
    }

    ;
    void printHelp() const throw() {
        cerr << "Program Options \n";
        cerr << "-help print these options\n";
        vector<Option>::const_iterator iter = options.begin();
        vector<Option>::const_iterator iend = options.end();
        while ( iter != iend ) {
            iter->writeToStream ( cerr );
            ++iter;
        }
        exit ( EXIT_FAILURE );
    }
    ;
};

}

#endif /* PROGRAMOPTIONS_H_ */
