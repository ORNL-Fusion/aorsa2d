
OBJ_DIR = ../obj
MOD_DIR = ../mod
EXE_DIR = exe
SRC_TST = src
SRC_DIR = ../src

F90 = gfortran
MOD_LOC = -J${MOD_DIR}
WARN = 
BOUNDS = -fbounds-check
DEBUG = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid
DOUBLE = -fdefault-real-8
F90FLAGS = ${MOD_LOC} ${BOUNDS} ${DEBUG} ${WARN}

NETCDF_DIR = /home/dg6/code/netcdf/netcdf_gnu64
NETCDF = -I ${NETCDF_DIR}/include -L ${NETCDF_DIR}/lib -lnetcdf

CPP_DIRECTIVES = -DcylXYZ -DzFunHammett

.PHONY : clean

all: ${EXE_DIR}/test_expBesselI \
		${EXE_DIR}/test_sigma

${EXE_DIR}/%: ${SRC_TST}/%.f90 
	${F90} ${F90FLAGS} $< -o $@ ${OBJ_DIR}/*.o ${NETCDF}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF} ${CPP_DIRECTIVES}

# Double

include ../Makefile.double

# Dependencies

${EXE_DIR}/test_expBesselI: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/bessel.o \
	${OBJ_DIR}/check_mod.o

${EXE_DIR}/test_sigma: \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/sigma.o \
		${OBJ_DIR}/check_mod.o

include ../Makefile.deps

clean:
	rm ${OBJ_DIR}/*.o ${MOD_DIR}/*.mod ${EXE_DIR}/* *.nc
