
OBJ_DIR = ../obj
MOD_DIR = ../mod
EXE_DIR = exe
SRC_TST = src
SRC_DIR = ../src

F90 = gfortran
MOD_LOC = -J${MOD_DIR}
WARN = 
BOUNDS = -fbounds-check
DEBUG = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid#,underflow,overflow
DOUBLE = -fdefault-real-8
F90FLAGS = ${MOD_LOC} ${BOUNDS} ${DEBUG} ${WARN}

NETCDF_DIR = ${HOME}/code/netcdf/gnu_4.3.2
NETCDF = -I ${NETCDF_DIR}/include -L ${NETCDF_DIR}/lib -lnetcdf -lnetcdff

CPP_DIRECTIVES = -DcylXYZ -DzFunHammett

.PHONY : clean

all: ${EXE_DIR}/test_expBesselI \
		${EXE_DIR}/test_sigma \
		${EXE_DIR}/test_chebyshev \
		${EXE_DIR}/test_cerf \
		${EXE_DIR}/test_spline

${EXE_DIR}/%: ${SRC_TST}/%.f90 
	${F90} ${F90FLAGS} $< -o $@ ${OBJ_DIR}/*.o ${NETCDF}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF} ${CPP_DIRECTIVES}

${OBJ_DIR}/%: ${SRC_TST}/cerf_wrapper.c 
	gcc -c $< -o $@ 


# Dependencies

${EXE_DIR}/test_expBesselI: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/bessel.o \
	${OBJ_DIR}/check_mod.o

${EXE_DIR}/test_sigma: \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/sigma.o \
		${OBJ_DIR}/check_mod.o

${EXE_DIR}/test_chebyshev: \
		${OBJ_DIR}/chebyshev_mod.o

${EXE_DIR}/test_cerf: \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/zfunHammett.o \
		${OBJ_DIR}/cerf_wrapper.o \
		${OBJ_DIR}/z_erfc.o

${EXE_DIR}/test_spline: \
		${OBJ_DIR}/spline_dlg.o \
		${OBJ_DIR}/check_mod.o

include ../Makefile.deps

# Double

include ../Makefile.double


clean:
	rm ${OBJ_DIR}/*.o ${MOD_DIR}/*.mod ${EXE_DIR}/* *.nc
