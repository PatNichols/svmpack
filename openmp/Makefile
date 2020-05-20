#CXX=g++
#CXXFLAGS=-O3 -funroll-all-loops -finline-limit=2000 -ffast-math -fstrict-aliasing -msseregparm -mfpmath=sse,387 \
#	-msse2 -msse3 -msse4a -m3dnow -m64 -minline-all-stringops -mtune=amdfam10  -I$(BOOST_DIR)/include -Wall
#CXXFLAGS=-ggdb -Wall -I/home/d3p708/include

CXXFLAGS= -O3 -I/Users/pnichols/homebrew/Cellar/gcc/9.3.0/lib/gcc/9/gcc/x86_64-apple-darwin18/9.3.0/include -ffast-math -march=native -mavx2 -fexpensive-optimizations -fopenmp

#CXXFLAGS= -O2 -funroll-all-loops -finline-limit=2000 -ffast-math \
#		-fstrict-aliasing -mfpmath=sse,387 -fno-rtti -m64 \
#	-msse2 -msse3 -msse4a -m3dnow -minline-all-stringops -mtune=core2  \
#	-fexpensive-optimizations -Wall 


#CXXFLAGS= -O3 -m64 -mtune=core2 -msse4a -ffast-math -mfpmath=sse,387 -march=core2
#	      -fprefetch-loop-arrays -funroll-all-loops -Wl,-z common-page-size=2M 

#CXXFLAGS=-ggdb -g 

LIBS= -L/Users/pnichols/homebrew/lib -lboost_thread-mt -lpthread -lm 
SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.h)
OBJS=$(patsubst %.cpp,%.o,$(SRCS))
    

all:  svmtask svmtranslate
        
svmtask: $(HDRS) $(OBJS) svmtask.cpp
	$(CXX) $(CXXFLAGS) -o svmtask svmtask.cpp $(LIBS)

svmtranslate: $(HDRS) $(OBJS) $(SRCS)
	$(CXX) $(CXXFLAGS) -o svmtranslate svmtranslate.cpp $(LIBS)

            
clean:
	rm -f svmtask
	rm -f svmtranslate
	rm -f *.o
                    
