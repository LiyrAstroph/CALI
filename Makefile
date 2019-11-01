all: cali

SYSTEM="Linux"
#SYSTEM="Darwin"

ifeq ($(SYSTEM), "Darwin")
FITS_INC = -I /opt/local/include/cfitsio
FITS_LIB = -L/opt/local/lib -lcfitsio
#PLPLOT_INC = -I /usr/local/include/plplot
#PLPLOT_LIB = -L /usr/local/lib -lplplotd -lm

GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I /usr/local/share/lapack/include -I/opt/local/include
LAPACK_LIBS = -framework vecLib -L /usr/local/share/lapack/lib -llapacke -llapack 
#-lcblas 
CBLAS_INCL  =
CBLAS_LIBS  =     
OPTIMIZE    = -O2
endif

ifeq ($(SYSTEM), "Linux")

GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl)

LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas -lgfortran

OPTIMIZE    = -O2 
endif

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL) $(FITS_INC)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS) $(PLPLOT_LIB) $(FITS_LIB)

OBJS = main.o allvars.o filerw.o calibration.o mathfun.o drw.o reconstruct.o

INCL = proto.h allvars.h version.h

$(OBJS): $(INCL)

cali: Makefile $(OBJS) $(INCL)
	gcc -o $@ $(OBJS) $(CFLAGS) $(LIBS)

clean:
	rm *.o cali
