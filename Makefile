all: cali

SYSTEM="Linux"
#SYSTEM="Darwin"

ifeq ($(SYSTEM), "Darwin")

GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I /usr/local/share/lapack/include -I/opt/local/include
LAPACK_LIBS = -framework vecLib -L /usr/local/share/lapack/lib -llapacke -llapack 
    
OPTIMIZE    = -O2
endif

ifeq ($(SYSTEM), "Linux")

GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl)

LAPACK_INCL = $(shell pkg-config --cflags lapacke lapack) 
LAPACK_LIBS = $(shell pkg-config --libs lapack lapacke)

OPTIMIZE    = -O2 
endif

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL) $(FITS_INC)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS) $(PLPLOT_LIB) $(FITS_LIB)

OBJS = main.o allvars.o filerw.o calibration.o mathfun.o drw.o reconstruct.o \
       run.o

INCL = proto.h allvars.h version.h

$(OBJS): $(INCL)

cali: Makefile $(OBJS) $(INCL)
	gcc -o $@ $(OBJS) $(CFLAGS) $(LIBS)

clean:
	rm *.o cali
