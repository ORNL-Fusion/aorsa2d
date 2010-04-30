
OBJ_DIR = obj
MOD_DIR = mod
EXE_DIR = exe
SRC_DIR = src
SRC_AORSA = ../src

F90 = gfortran
MOD_LOC = -J${MOD_DIR}
WARN = 
BOUNDS = -fbounds-check
DEBUG = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid
DOUBLE = -fdefault-real-8
F90FLAGS = ${MOD_LOC} ${BOUNDS} ${DEBUG} ${WARN}

NETCDF_DIR = /home/dg6/code/netcdf/netcdf_gnu64
NETCDF = -I ${NETCDF_DIR}/include -L ${NETCDF_DIR}/lib -lnetcdf

.PHONY : clean

${EXE_DIR}/test_expBesselI: ${SRC_DIR}/test_expBesselI.f90 
	${F90} ${F90FLAGS} $< -o $@ ${OBJ_DIR}/*.o ${NETCDF}

${OBJ_DIR}/%.o: ${SRC_AORSA}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF}

${OBJ_DIR}/bessel.o: ${SRC_AORSA}/bessel.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${DOUBLE}

# Dependencies

${EXE_DIR}/test_expBesselI: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/bessel.o \
	${OBJ_DIR}/check.o

clean:
	rm ${OBJ_DIR}/*.o ${MOD_DIR}/*.mod ${EXE_DIR}/* *.nc
