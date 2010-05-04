.SUFFIXES : .o .f90 .f .F90

EXEC = xaorsa2d
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

FFT_DIR = $(SRC_DIR)/fftpack
CQL_DIR = $(SRC_DIR)/cql3d
INCLUDE_DIR = $(SRC_DIR)/cql3d
DLG_DIR = $(SRC_DIR)/dlg


# objects
# -------

OBJ_FILES := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*.*)))
OBJ_DLG := $(patsubst src/dlg/%,obj/%.o,$(basename $(wildcard src/dlg/*)))

#OBJ_CQL3D := $(patsubst src/cql3d/%,obj/%.o,$(basename $(wildcard src/cql3d/*.*)))
#OBJ_FFT := $(patsubst src/fftpack/%,obj/%.o,$(basename $(wildcard src/fftpack/*)))

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

# set the MODE to "serial" or "parallel"

MODE = "parallel"

# set solve precision to "single" or "double" 

SOLVE_PREC = "double"

# set the coordinate system to an understandable 
# cylindrical ("cylProper") or old version ("cylXYZ")

COORDSYS = "cylXYZ"


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


# compile flags
# -------------

BOUNDS = -fbounds-check
WARN = #-Wall
DEBUG = -pg -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid#,overflow#,underflow
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
	LIBS = ${BLAS} ${LAPACK} ${NETCDF} ${PAPI}
endif
INC_DIR = ${PAPI_INC}

F90FLAGS = ${WARN} ${DEBUG} ${BOUNDS} ${MOD_LOC} 
LD_FLAGS =

.PHONY: depend clean

$(EXEC): $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) 
	$(F90) -o $(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D) $(OBJ_DLG) $(LIBS) ${LINK_FLAGS}
			    		    		     		   			     			     				
## FFT files
#
#${OBJ_DIR}/%.o: ${FFT_DIR}/%.f
#	${F90} -c ${F90FLAGS} $< -o $@ 
				
# SRC files

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR}

## CQL files
#
#${OBJ_DIR}/%.o: ${CQL_DIR}/%.*
#	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF}

# DLG files
		
${OBJ_DIR}/%.o: ${DLG_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF} ${INC_DIR} ${CPP_DIRECTIVES}

${OBJ_DIR}/%.o: ${DLG_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF} ${INC_DIR} ${CPP_DIRECTIVES}


# Double precision routines
# -------------------------

${OBJ_DIR}/bessel.o: ${SRC_DIR}/bessel.f90
	${F90} -c ${F90FLAGS} ${DOUBLE} $< -o $@ ${NETCDF} ${INC_DIR}

${OBJ_DIR}/zfunction.o: ${SRC_DIR}/zfunction.f90
	${F90} -c ${F90FLAGS} ${DOUBLE} $< -o $@ ${NETCDF} ${INC_DIR}

${OBJ_DIR}/Zfun.o: ${SRC_DIR}/Zfun.f90
	${F90} -c ${F90FLAGS} ${DOUBLE} $< -o $@ ${NETCDF} ${INC_DIR}

${OBJ_DIR}/ztable.o: ${SRC_DIR}/ztable.f90
	${F90} -c ${F90FLAGS} ${DOUBLE} $< -o $@ ${NETCDF} ${INC_DIR}

# Dependencies	

${OBJ_DIR}/aorsa2dMain.o: \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/eqdsk_dlg.o \
		${OBJ_DIR}/sigma.o \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/interp.o \
		${OBJ_DIR}/inv_fourier.o \
		${OBJ_DIR}/write_data.o \
		${OBJ_DIR}/grid.o \
		${OBJ_DIR}/bField.o \
		${OBJ_DIR}/profiles.o \
		${OBJ_DIR}/rotation.o \
		${OBJ_DIR}/mat_fill.o \
		${OBJ_DIR}/antenna.o \
		${OBJ_DIR}/solve.o \
		${OBJ_DIR}/timer.o

${OBJ_DIR}/eqdsk_dlg.o: \
		${OBJ_DIR}/dlg.o \
		${OBJ_DIR}/fitpack.o \
		${OBJ_DIR}/grid.o \
		${OBJ_DIR}/aorsa2din_mod.o

${OBJ_DIR}/sigma.o: \
		${OBJ_DIR}/bessel.o \
		${OBJ_DIR}/zfunction.o \
		${OBJ_DIR}/bField.o \
		${OBJ_DIR}/rotation.o

${OBJ_DIR}/zfunction.o: \
		${OBJ_DIR}/ztable.o \
		${OBJ_DIR}/Zfun.o

${OBJ_DIR}/interp.o: \
		${OBJ_DIR}/fitpack.o \
		${OBJ_DIR}/eqdsk_dlg.o

${OBJ_DIR}/Zfun.o: \
		${OBJ_DIR}/constants.o

${OBJ_DIR}/bField.o: \
		${OBJ_DIR}/interp.o \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/grid.o

${OBJ_DIR}/profiles.o: \
		${OBJ_DIR}/bField.o \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/aorsa2din_mod.o

${OBJ_DIR}/rotation.o: \
		${OBJ_DIR}/bField.o \
		${OBJ_DIR}/derivatives.o \
		${OBJ_DIR}/grid.o \
		${OBJ_DIR}/eqdsk_dlg.o

${OBJ_DIR}/grid.o: \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/parallel.o

${OBJ_DIR}/mat_fill.o: \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/sigma.o \
		${OBJ_DIR}/grid.o \
		${OBJ_DIR}/rotation.o \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/profiles.o \
		${OBJ_DIR}/bField.o \
		${OBJ_DIR}/parallel.o \
		${OBJ_DIR}/eqdsk_dlg.o

${OBJ_DIR}/antenna.o: \
		${OBJ_DIR}/grid.o \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/profiles.o \
		${OBJ_DIR}/parallel.o \
		${OBJ_DIR}/eqdsk_dlg.o

${OBJ_DIR}/write_data.o: \
		${OBJ_DIR}/mat_fill.o \
		${OBJ_DIR}/solve.o \
		${OBJ_DIR}/constants.o

${OBJ_DIR}/inv_fourier.o: \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/grid.o

${OBJ_DIR}/solve.o: \
		${OBJ_DIR}/parallel.o

${OBJ_DIR}/write_data.o: \
		${OBJ_DIR}/check_mod.o

clean:
	rm $(EXEC) $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod


