CXX=g++
#CXXFLAGS=-O3 -funroll-all-loops -finline-limit=2000 -ffast-math -fstrict-aliasing -msseregparm -mfpmath=sse,387 \
#	-msse2 -msse3 -msse4a -m3dnow -m64 -minline-all-stringops -mtune=amdfam10  -I$(BOOST_DIR)/include -Wall
#CXXFLAGS=-ggdb -Wall -I/home/d3p708/include

CXXFLAGS= -O3 -ffast-math \
		-fstrict-aliasing -mfpmath=sse -m64 \
	-mavx -minline-all-stringops -mtune=corei7-avx \
	-fexpensive-optimizations 

#CXXFLAGS= -O2 -funroll-all-loops -finline-limit=2000 -ffast-math \
#		-fstrict-aliasing -mfpmath=sse,387 -fno-rtti -m64 \
#	-msse2 -msse3 -msse4a -m3dnow -minline-all-stringops -mtune=core2  \
#	-fexpensive-optimizations -Wall 


#CXXFLAGS= -O3 -m64 -mtune=core2 -msse4a -ffast-math -mfpmath=sse,387 -march=core2
#	      -fprefetch-loop-arrays -funroll-all-loops -Wl,-z common-page-size=2M 

#CXXFLAGS=-ggdb -g 

LIBS= -L/usr/lib/x86_64_linux_gnu -lboost_thread -lpthread -lm -lrt
SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.h)
OBJS=$(patsubst %.cpp,%.o,$(SRCS))
    
%.o:%.cpp  $(HDRS) $(SRCS)
	$(CXX) -c $(CXXFLAGS) $< -o $@
        
all:  svmtask
        
svmtask: $(HDRS) $(OBJS) $(SRCS)
	$(CXX) $(CXXFLAGS) -o svmtask $(OBJS) $(LIBS)
            
clean:
	rm -f svmtask
	rm -f *.o
                    
