GSL_DIR=/software/groups/ping_group/shared/libs-CentOS9/gsl-2.7.1/build
MKLROOT=/software/groups/ping_group/shared/libs-CentOS9/mkl-2024.1.0/mkl/2024.1
IFLAGS=-I${GSL_DIR}/include -I${MKLROOT}/include
CC=mpicxx -std=c++11 -O2 -g 
GSL_LIBS=-L${GSL_DIR}/lib -lgsl -lgslcblas
LAPACK_LIBS=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LDLIBS=${GSL_LIBS} ${LAPACK_LIBS}
LDFLAGS = -Wl,-rpath,${GSL_DIR}/lib,-rpath,${MKLROOT}/lib/intel64

