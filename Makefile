.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# .PHONY : usage
# usage:
# 	@echo ""
# 	@echo "    * USAGE * "
# 	@echo ""
# 	@echo ""

objdir = .obj

# Command-line options at make call
ifort ?= 0
debug ?= 1

FC = gfortran

LIBSYS = /usr/lib
INCSYS = /usr/include

LIBNETCDF = -I$(INCSYS) -L$(LIBSYS) -lnetcdff -lnetcdf
LIBNCIO   = -Wl,-rpath=../ncio -I../ncio/.obj -L../ncio -lncio

SRCFAGS =
OBJFAGS = -J$(objdir) -I$(objdir)
BINFAGS =

# LIBNCIO   = -Wl,rpath=../NCIO/ncio/ -L../NCIO/ncio -lncio
# LIBNCIO = -I/home/dmr/Software-Devel/NCIO/ncio/.obj # -L ../NCIO/ncio/ -lncio

LIBS= $(LIBNCIO) $(LIBNETCDF)

DFLAGS = -O3 -cpp -ffixed-line-length-132 -fno-align-commons -fdefault-real-8 -fdefault-double-8
ifeq ($(debug), 1)
    DFLAGS  = -g -pg -cpp -Wall -ffixed-line-length-132 -fno-align-commons -fdefault-real-8 -fdefault-double-8 -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
endif

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
    FC = ifort
    LIBNETCDF = -I /usr/local/install/netcdf-4.1.1/include -L/usr/local/install/netcdf-4.1.1/lib -lnetcdf -lnetcdff  -L/lib64

	FLAGS        = -module $(objdir) -L$(objdir) -I$(INC)
	DFLAGS       = -cpp -vec-report0 -132
	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	    # -w
	endif
    OBJFAGS = -I$(objdir)
endif
# 


SRC     = \
        Grid_utils.f90 \
        Grid_Class.f90 \
        Sub_Grid_Class.f90\
        toms760-wrapper.f90 \
        distance_great_circle.f90\
        interpolate_mod.f90\
        module_test.f90

OBJS    = $(SRC:.f90=.o)
MODS    = ${SRC:.f90=.mod}
OBJPATH = $(addprefix $(objdir)/, $(OBJS))

## Individual libraries or modules ##
# %.mod %.o : %.f90
%.o : %.f90
	$(FC) $(DFLAGS) $(FLAGS) ${LIBS} $(OBJFAGS) -o $(objdir)/$@ -c $<
	@echo "" $(OBJPATH)

%.mod : %.f90 %.o
	@true

## Shared library 
# .so: %.f90
# 	$(FC) -c -shared -fPIC $(DFLAGS) $(FLAGS) ${LIBS} -o ncio.so $<

## Complete programs

main: $(OBJS)
	$(FC) $(DFLAGS) -o test.x $(OBJFAGS) $(OBJPATH) test.f90 $(LIBS)
	@echo " "
	@echo "    test.x is ready."
	@echo " "
lib:  $(OBJS)
	ar -cvr libgrid.a $(OBJPATH)
	ranlib libgrid.a
	mv libgrid.a library/.
	@echo " "
	@echo "    libgrid.a is ready."
	@echo " "
clean:
	rm -f test.x $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.so gmon.out
