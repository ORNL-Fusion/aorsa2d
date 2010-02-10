HOME = /home/dg6/code
EXEC = xaorsa2d
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

FFT_DIR = $(SRC_DIR)/fftpack
CQL_DIR = $(SRC_DIR)/cql3d
INCLUDE_DIR = $(SRC_DIR)/cql3d
DLG_DIR = $(SRC_DIR)/dlg

OBJ_FFT := $(patsubst src/fftpack/%,obj/%.o,$(basename $(wildcard src/fftpack/*)))
OBJ_DLG := $(patsubst src/dlg/%,obj/%.o,$(basename $(wildcard src/dlg/*)))
OBJ_FILES := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*.*)))
OBJ_CQL3D := $(patsubst src/cql3d/%,obj/%.o,$(basename $(wildcard src/cql3d/*.*)))

BLACS = 
BLAS = 
SCALAPACK = 
NETCDF_DIR = /home/dg6/code/netcdf/netcdf_gnu64
NETCDF = -I ${NETCDF_DIR}/include -L ${NETCDF_DIR}/lib -lnetcdf 

LIBS = $(SCALAPACK) $(BLACS) $(BLAS) ${NETCDF}
INC_DIR = 

BOUNDS = -fbounds-check
WARN = #-Wall
DEBUG = -pg -g -fbacktrace -fsignaling-nans #-ffpe-trap=zero,invalid#,overflow#,underflow
F90 = gfortran
LINK = gfortran
MOD_LOC = -Jmod
F90FLAGS = ${WARN} ${DEBUG} ${BOUNDS} ${MOD_LOC}
LD_FLAGS =

.PHONY: depend clean

$(EXEC): $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) 
	$(LINK) -o $(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D) $(OBJ_DLG) $(LIBS) ${LINK_FLAGS}
			    		    		     		   			     			     				
# FFT files

${OBJ_DIR}/%.o: ${FFT_DIR}/%.*
	${F90} -c ${F90FLAGS} $< -o $@ 
				
# SRC files

${OBJ_DIR}/%.o: ${SRC_DIR}/%.*
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR}

# CQL files

${OBJ_DIR}/%.o: ${CQL_DIR}/%.*
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF}

# DLG files
		
${OBJ_DIR}/%.o: ${DLG_DIR}/%.*
	${F90} -c ${F90FLAGS} $< -o $@ ${NETCDF} ${INC_DIR} ${CPP_DIRECTIVES}

# Dependencies	

${OBJ_DIR}/aorsa2dMain.o: \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/eqdsk_dlg.o \
		${OBJ_DIR}/aorsaSubs.o \
		${OBJ_DIR}/sigma.o \
		${OBJ_DIR}/aorsa2din_mod.o \
		${OBJ_DIR}/interp.o

${OBJ_DIR}/eqdsk_dlg.o: \
		${OBJ_DIR}/dlg.o \
		${OBJ_DIR}/fitpack.o

${OBJ_DIR}/aorsaSubs.o: \
		${OBJ_DIR}/bessel.o

${OBJ_DIR}/sigma.o: \
		${OBJ_DIR}/bessel.o \
		${OBJ_DIR}/zfunction.o

${OBJ_DIR}/zfunction.o: \
		${OBJ_DIR}/ztable.o \
		${OBJ_DIR}/aorsaSubs.o

${OBJ_DIR}/interp.o: \
		${OBJ_DIR}/fitpack.o \
		${OBJ_DIR}/eqdsk_dlg.o

clean:
	rm $(EXEC) $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod


