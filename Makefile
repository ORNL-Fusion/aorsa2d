.SUFFIXES : .o .f90 .f .F90

AORSA2D = xaorsa2d
AORSA1D = xaorsa1d
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

# objects
# -------

#OBJ_FILES := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*.*)))
OBJ_FILES = $(wildcard obj/*.o)

# libraries 
# ---------

BLAS = ${HOME}/code/goto_blas/GotoBLAS2/libgoto2.a -pthread
LAPACK = ${HOME}/code/lapack/lapack-3.1.1/lapack_LINUX.a 
NETCDF_DIR = /home/dg6/code/netcdf/netcdf_gnu64
NETCDF = -I ${NETCDF_DIR}/include -L ${NETCDF_DIR}/lib -lnetcdf 
BLACS = \
	${HOME}/code/blacs/blacs_gnu64/LIB/blacs_MPI-LINUX-0.a \
	${HOME}/code/blacs/blacs_gnu64/LIB/blacsF77init_MPI-LINUX-0.a \
	${HOME}/code/blacs/blacs_gnu64/LIB/blacs_MPI-LINUX-0.a
SCALAPACK = ${HOME}/code/scalapack/scalapack_gnu64/libscalapack.a
PAPI_INC = -I/usr/include 
PAPI = -lpapi
CUDA_DIR = ${HOME}/code/cuda/3.0/cuda
CUDA = -L ${CUDA_DIR}/lib64 -lcublas -lcudart -lcuda -L /usr/lib64
MAGMA_DIR = ${HOME}/code/magma/magma_0.2
MAGMA = -L ${MAGMA_DIR}/lib -lmagma -lmagmablas ${MAGMA_DIR}/lib/libmagma_64.a

# set the MODE to "serial" or "parallel"

MODE = "serial"

# set solve precision to "single" or "double" 

SOLVE_PREC = "double"

# set the coordinate system to an understandable 
# cylindrical ("cylProper") or old version ("cylXYZ")

COORDSYS = "cylXYZ"

# set the z function to use
# either "zFunOriginal" or "zFunHammett"

ZFUN = "zFunHammett"

# use the Chebychev basis set

BASIS = "chebychev"


# pre-processor directives
# ------------------------

ifeq (${MODE},"parallel")
	CPP_DIRECTIVES = -Dpar
endif

ifeq (${SOLVE_PREC},"double")
	CPP_DIRECTIVES := -Ddblprec ${CPP_DIRECTIVES}
endif

ifeq (${COORDSYS},"cylProper")
	CPP_DIRECTIVES := -DcylProper ${CPP_DIRECTIVES}
endif

ifeq (${ZFUN},"zFunHammett")
	CPP_DIRECTIVES := -DzFunHammett ${CPP_DIRECTIVES}
endif

ifeq (${BASIS},"chebychev")
	CPP_DIRECTIVES := -Dchebychev ${CPP_DIRECTIVES}
endif

# compile flags
# -------------

BOUNDS = -fbounds-check 
WARN = #-Wall
DEBUG = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid#,overflow#,underflow
DOUBLE = -fdefault-real-8
ifeq (${MODE},"parallel")
	F90 = mpif90
else
	F90 = gfortran
endif
MOD_LOC = -Jmod


# other machines
# --------------

ifeq (${HOME},/global/homes/g/greendl1)
	include Makefile.franklin
endif

ifeq (${HOME},/ccs/home/dg6)
	include Makefile.jaguarpf
endif

# the order of linking libs is important
ifeq (${MODE},"parallel")
	LIBS = ${SCALAPACK} ${BLACS} ${BLAS} ${LAPACK} ${NETCDF} ${PAPI}
else
	LIBS = ${LAPACK} ${NETCDF} ${PAPI} ${MAGMA} ${CUDA} -lstdc++ ${BLAS} 
endif
INC_DIR = ${PAPI_INC}

F90FLAGS = ${WARN} ${DEBUG} ${BOUNDS} ${MOD_LOC} 
LINK_FLAGS = 

.PHONY: depend clean

${AORSA2D}: ${SRC_DIR}/aorsa2dMain.F90  
	$(F90) ${F90FLAGS} ${SRC_DIR}/aorsa2dMain.F90 -o ${AORSA2D} $(OBJ_FILES) $(LIBS) ${LINK_FLAGS} ${CPP_DIRECTIVES} ${INC_DIR}

${AORSA1D}: ${SRC_DIR}/aorsa1dMain.F90  
	$(F90) ${F90FLAGS} ${SRC_DIR}/aorsa1dMain.F90 -o ${AORSA1D} $(OBJ_FILES) $(LIBS) ${LINK_FLAGS} ${CPP_DIRECTIVES} ${INC_DIR}
			    		    		     		   			     			     				
# SRC files

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR}

#${OBJ_DIR}/solve.o: ${SRC_DIR}/solve.F90
#	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR} -fno-underscoring
#
#${OBJ_DIR}/get_nb.o: ${MAGMA_DIR}/testing/get_nb.cpp
#	${F90} -c $< -o $@
#
#${OBJ_DIR}/fortran.o: ${CUDA_DIR}/src/fortran.c
#	gcc -c $< -o $@ -I ${CUDA_DIR}/include


# Double precision routines
# -------------------------

include Makefile.double

# Dependancies
# ------------

include Makefile.deps

clean:
	rm ${AORSA1D} ${AORSA2D} $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod


