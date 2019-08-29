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
#FITS_INC = -I /usr/include/cfitsio
#FITS_LIB = -lcfitsio
#PLPLOT_INC = -I/usr/include/plplot
#PLPLOT_LIB = -lplplotd -lm

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

filerw.o: filerw.c allvars.h Makefile
	gcc -c -o $@ $< $(CFLAGS) 

allvars.o: allvars.c allvars.h Makefile
	gcc -c -o $@ $<	$(CFLAGS) 
	
calibration.o: calibration.c allvars.c Makefile
	gcc -c -o $@ $< $(CFLAGS)

mathfun.o: mathfun.c allvars.h Makefile
	gcc -c -o $@ $< $(CFLAGS)

drw.o: drw.c allvars.h Makefile
	gcc -c -o $@ $< $(CFLAGS)

reconstruct.o: reconstruct.c allvars.h Makefile
	gcc -c -o $@ $< $(CFLAGS)

main.o: main.c allvars.h Makefile
	gcc -c -o $@ $< $(CFLAGS)

#plot: clea_plot.c Makefile 
#	gcc -o $@ $< $(CFLAGS) $(LIBS)

cali: Makefile $(OBJS) proto.h
	gcc -o $@ $(OBJS) $(CFLAGS) $(LIBS)

clean:
	rm *.o cali
