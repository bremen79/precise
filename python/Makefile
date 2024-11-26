CXX = g++
CXX_OPTIONS = -std=c++0x -O3 -msse -msse2 -ffast-math -DNDEBUG
INCLUDE = -I.

opt.o: opt.cpp
	$(CXX) $(INCLUDE) $(CXX_OPTIONS) -o opt.o -c -fpic opt.cpp

opt_ctypes: opt.o
	$(CXX) $(INCLUDE) $(CXX_OPTIONS) opt.o -o lib_opt.so --shared

all: opt.o opt_ctypes

clean:
	rm *.o; rm *.so; rm *.pyc
