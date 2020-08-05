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
#include <cstring>
#include <unistd.h>
#include "parse.h"
namespace svmpack
{


class ProgramOptions
{
    struct Option {
        explicit Option ( const char *keyword ) :
                key ( keyword ), des (), val ( string ( "" ) ), stat ( 0 ) {
        }
        

        explicit Option ( const char *keyword, const char *description ) :
                key ( keyword ), des ( description ), val ( string ( "" ) ), stat ( 0 ) {
        }
        

        explicit Option ( const char *keyword, const char *description, const char *value ) :
                key ( keyword ), des ( description ), val ( value ), stat ( -1 ) {
        }
        
        explicit Option ( const std::string& keyword ) :
                key ( keyword ), des (), val ( std::string ( "" ) ), stat ( 0 ) {
        }
        

        explicit Option ( const std::string& keyword, const std::string& description ) :
                key ( keyword ), des ( description ), val ( std::string ( "" ) ), stat ( 0 ) {
        }
        

        explicit Option ( const std::string& keyword, const std::string& description, const std::string& value ) :
                key ( keyword ), des ( description ), val ( value ), stat ( -1 ) {
        }
        

        bool hasValue() const throw() {
            return ( stat != 0 );
        }

        bool matches ( const char *keyword ) const throw() {
            return ( key.compare ( keyword ) == 0 );
        }
        

        string getValue() const throw () {
            if ( stat ) {
                return string ( val );
            }
            cerr << " value requested but none was set in ProgramOption::Option::getValue()\n";
            exit ( EXIT_FAILURE );
        }
        

        void setValue ( const char *valueString ) throw () {
            if ( stat != 1 ) {
                val = valueString;
            }
            stat = 1;
        }
        

        ostream& writeToStream ( ostream& os ) const throw() {
            os << "-" << key << " " << des << "\n";
            if ( stat == -1 ) {
                os << "\t\t\t\t (default = " << val << ")\n";
            }
            return os;
        }

        string key;
        string des;
        string val;
        int stat;
    };

    vector<Option> options;
    std::string prog_name;
public:
    ProgramOptions() {
    }
    ;

    ~ProgramOptions() {
    }
    ;

    void toUpper(std::string& s) {
        for (std::size_t k=0;k<s.size();++k) {
            s[k] = toupper(s[k]);
        }
    }

    void parseEnv ( const char * Prefix = 0x0 ) noexcept {
        vector<Option>::iterator iter = options.begin();
        vector<Option>::iterator iend = options.end();
        if (!Prefix) {
          while (iter!=iend) {
            std::string key = iter -> key;
            toUpper(key);
            char * p = getenv(key.c_str());
            if (p) {
                setValue( (iter->key).c_str(), p);
            }
            ++iter;
          }
          return;        
        }
        while (iter!=iend) {
            std::string s_key = iter -> key;
            toUpper(s_key);
            std::string key = std::string(Prefix) + "_" + s_key;
            char * p = getenv(key.c_str());
            if (p) {
                setValue( (iter->key).c_str(), p);
            }
            ++iter;
        }    
    }


    void parseInputFile  ( const char *filename ) throw() {
        std::ifstream in(filename);
        std::string sline;
        std::string delims(" :=\n");
        std::vector<std::string> tokens;
        if ( !in ) {
            cerr << "could not open the input options file " << filename << "\n";
            exit ( EXIT_FAILURE ); 
        }
        while ( in ) {
            getline(in,sline);
            if (sline.size() && sline[0]=='#') continue;
            explodeString(sline,delims,tokens);          
            auto nsz = tokens.size();
            if (nsz==0) continue;  // blank line allowed
            if (nsz==1) {
                options.push_back(Option(tokens[0]));
            } else {
                if (nsz==2) {
                    options.push_back(Option(tokens[0],tokens[1]));
                }else{
                    if (nsz==3) {
                        options.push_back(Option(tokens[0],tokens[1],tokens[2]));
                    }else{
                        std::cerr << "Bad format in options input file\n";
                        exit(EXIT_FAILURE);
                    }
                }
            } 
            if (in.eof()) break;
        }
        in.close();    
    }

    void parseConfigFile ( const char *filename ) throw() {
        std::ifstream in(filename);
        std::string sline;
        std::string delims(" :=\n");
        std::vector<std::string> tokens;
        if ( !in ) {
            cerr << "could not open the options config file " << filename << "\n";
            exit ( EXIT_FAILURE ); 
        }
        while ( in ) {
            getline(in,sline);
            if (sline.size() && sline[0]=='#') continue; // comment lines allowed
            explodeString(sline,delims,tokens);          
            auto nsz = tokens.size();
            if (nsz==0) continue;  // blank line allowed
            if (nsz==1) {
                setValue(tokens[0].c_str(),"true");
            } else {
                if (nsz==2) {
                    setValue(tokens[0].c_str(),tokens[1].c_str());
                }else{
                    std::cerr << "Bad format in options config file\n";
                    exit(EXIT_FAILURE);
                }
            } 
            if (in.eof()) break;
        }
        in.close();    
     }


    void parseCommandLine ( int argc, char **argv ) throw() {
        prog_name = argv[0];
        for ( int k = 1; k < argc; ++k ) {
            std::string keystr = argv[k];
            if (keystr.size()>2) {
                std::size_t p = 1;
                if (keystr[0]!='-') {
                    std::cerr << "option expected but found " << keystr << "\n"; 
                    printHelp();
                }
                if (keystr[1]=='-') p = 2;
                std::string key = keystr.substr(p,std::string::npos);
                p = key.find_first_of("=",0);
                if (p==std::string::npos) {
                    // no equal sign found in argv
                    if (k!=(argc-1)) {
                        std::string val(argv[k+1]);
                        if (val.find_first_of("-")==std::string::npos) {
                            ++k;
                            setValue(key.c_str(),argv[k]);
                        } else{
                            // next word on command line is another option so we assume that the 
                            // option is implicitly true    
                            setValue(key.c_str(),"true");
                        }
                    }else{    
                            // next word on command line is followed by nothing so we assume that the 
                            // option is implicitly true    
                            setValue(key.c_str(),"true");                        
                    }            
                }else{
                    // we have -key=value or --key=value
                    std::string key_ = key.substr(0,p);
                    std::string val_ = key.substr(p+1,std::string::npos);
                    setValue(key_.c_str(),val_.c_str());
                }           
            }else{
                std::cerr << "Expected an option but found " << keystr << "\n";
                printHelp();
            } 
        }
    }


    void addOption ( const char *keyword, const char *description ) {
        options.push_back ( Option ( keyword, description ) );
    }


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


    template<class T> inline T getValue ( const char *keyword ) const throw() {
        vector<Option>::const_iterator iter = options.begin();
        vector<Option>::const_iterator iend = options.end();
        while ( iter != iend ) {
            if ( iter->matches ( keyword ) ) {
                return parse<T> ( iter->getValue() );
            }
            ++iter;
        }
        cerr << "could not find the option " << keyword << endl;
        exit ( EXIT_FAILURE );
    }


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

    void printHelp() const throw() {
        cerr << "Options for "<< prog_name << " are:\n";
        cerr << "-help print these options\n";
        vector<Option>::const_iterator iter = options.begin();
        vector<Option>::const_iterator iend = options.end();
        while ( iter != iend ) {
            iter->writeToStream ( cerr );
            ++iter;
        }
        exit ( EXIT_FAILURE );
    }

};

}

#endif /* PROGRAMOPTIONS_H_ */
