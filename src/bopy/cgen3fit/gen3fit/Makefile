CXX=g++
CXXFLAGS=-g -O2 -fPIC -std=c++0x -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5 -ffast-math
INCLUDE=-I/usr/include/python2.6 -I/usr/local/lib/python2.6/dist-packages/numpy/core/include
LIBS=-lgsl -lgslcblas

gen3fit.so: gen3fit.o cgen3fit.o clinfit.o
    ${CXX} ${CXXFLAGS} -shared -o gen3fit.so *.o ${LIBS}

gen3fit.o: gen3fit.cc cgen3fit.h cgen3fit.o clinfit.h clinfit.o
    ${CXX} ${CXXFLAGS} gen3fit.cc -c ${INCLUDE} ${LIBS}

cgen3fit.o: cgen3fit.cc cgen3fit.h
    ${CXX} ${CXXFLAGS} cgen3fit.cc -c ${LIBS}

clinfit.o: clinfit.cc clinfit.h
    ${CXX} ${CXXFLAGS} clinfit.cc -c ${LIBS}

clean: 
    rm -rf *.o *.so

