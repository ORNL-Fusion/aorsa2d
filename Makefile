.SUFFIXES : .o .f90 .f .F90

AORSA2D = xaorsa2d
AORSA1D = xaorsa1d
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod
CPP_DIR = cpp

COMPILER := GNU# GNU, PGI
PARALLEL := 0# 0, 1
GPU := 0# 0, 1
DDT := 0# 0, 1

ifeq (${GPU},1)
	COMPILER := GNU
endif

# objects
# -------

OBJ_FILES = $(wildcard obj/*.o)

# libraries 
# ---------
#
LIBS := 
INC_DIR := 
CPP_DIRECTIVES :=

# set solve precision to double 
# -----------------------------
CPP_DIRECTIVES += -Ddblprec 

# set the coordinate system to an understandable 
# cylindrical ("-DcylProper")  
# ----------------------------------------------
CPP_DIRECTIVES += -DcylProper

# set the z function to use
# either "zFunOriginal" or "zFunHammett"
# --------------------------------------
CPP_DIRECTIVES += -DzFunHammett 

# use papi, set to "-Dusepapi"
# ----------------------------
CPP_DIRECTIVES += #-Dusepapi 

# caculate sigma as part of the fill (=2) or standalone with file write (=1)
# --------------------------------------------------------------------------
CPP_DIRECTIVES += -D__sigma__=2

# Try removing the U rotation matricies (actually just set to the identity in places)
# and just have the output of the sigma routines be in rtz coords and not in alp,bet,prl.
# This is in an attempt to see how this is related to the noise problem with the poloidal
# field. 
# ---------------------------------------------------------------------------------------
CPP_DIRECTIVES += -D__noU__=0

# debug flags
# -----------

CPP_DIRECTIVES += -D__CheckParallelLocation__=0# Double check that myRow==pr_sp .and. myCol==pc_sp
CPP_DIRECTIVES += -D__debugSigma__=0
CPP_DIRECTIVES += -D__DebugSetMetal__=0
CPP_DIRECTIVES += -D__DebugBField__=0
CPP_DIRECTIVES += -D_DEBUG_ROTATION=0
CPP_DIRECTIVES += -DPRINT_SIGMA_ABP=0
CPP_DIRECTIVES += -DPRINT_SIGMA=0

# compile flags
# -------------

ifeq (${COMPILER},PGI)
    MOD_LOC:= -module mod
    DOUBLE:= -Mr8
    WARN:=
    DEBUG:= #-g -traceback -Ktrap=divz,inv,ovf
    OPTIMIZATION:= -fast
    BOUNDS:= #-Mbounds
    FORMAT:=
else
	FORMAT := -ffree-line-length-none
	BOUNDS := -fbounds-check 
	WARN := -Wall
	DEBUG := -pg -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid#,overflow#,underflow
	OPTIMIZATION := #-O3
	DOUBLE := -fdefault-real-8
	MOD_LOC := -Jmod
endif

ifeq (${PARALLEL},1)
    CPP_DIRECTIVES += -Dpar
endif


# Machine specific stuff
# ----------------------
ThisMachine := $(shell uname -n)

ifneq (,$(findstring titan,$(ThisMachine)))
ThisMachine := titan
endif
ifneq (,$(findstring hopper,$(ThisMachine)))
ThisMachine := hopper
endif
ifneq (,$(findstring edison,$(ThisMachine)))
ThisMachine := edison
endif

include machine_makefiles/Makefile.${ThisMachine}

ifeq (${DDT},1)
    LIBS+=-Bddt
endif

F90FLAGS = ${FORMAT} ${WARN} ${DEBUG} ${BOUNDS} ${MOD_LOC} ${OPTIMIZATION}
LINK_FLAGS = 

.PHONY: depend clean

${AORSA2D}: ${SRC_DIR}/aorsa.F90  
	$(F90) ${F90FLAGS} ${SRC_DIR}/aorsa.F90 -o \
			${AORSA2D} \
			$(OBJ_FILES) \
			${GPTL} $(LIBS) \
			${LINK_FLAGS} \
			${CPP_DIRECTIVES} \
			${INC_DIR}

# SRC files
#
${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${CPP_DIRECTIVES} ${INC_DIR}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${CPP_DIRECTIVES} ${INC_DIR}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${CPP_DIRECTIVES} ${INC_DIR}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${CPP_DIRECTIVES} ${INC_DIR}

${OBJ_DIR}/z_erf.o: ${SRC_DIR}/z_erf.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${INC_DIR} ${DOUBLE}

# Uncomment the following 4 rules for MAGMA implementation
# in addition to two lines in Makefile.deps and adding two 
# trailing underscores in src/solve.f90

#${OBJ_DIR}/solve.o: ${SRC_DIR}/solve.F90
#	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF} ${CPP_DIRECTIVES} ${INC_DIR} -fno-underscoring
#
#${OBJ_DIR}/get_nb.o: ${MAGMA_DIR}/testing/get_nb.cpp
#	${F90} -c $< -o $@ -fPIC
#
#${OBJ_DIR}/fortran.o: ${CUDA_DIR}/src/fortran.c
#	gcc -c $< -o $@ -I ${CUDA_DIR}/include
#
#${OBJ_DIR}/magma_solve.o: ${SRC_DIR}/magma_solve.cpp 
#	gcc -c -g $< -o $@ -I ${CUDA_DIR}/include -I ${MAGMA_DIR}/include


# Double precision routines
# -------------------------

include Makefile.double

# Dependancies
# ------------

include Makefile.deps

clean:
	rm -f ${AORSA1D} ${AORSA2D} $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod aorsa.o *.mod


