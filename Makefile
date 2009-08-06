#
#	Makefile for Jaguar xt4, Franklin xt4 & dlghp
#
#	IMPORTANT: load modules with 
#		"module load netcdf" 
#		"module load pgplot"
#
#	For performance data, load the following before compiling:
#		"module load papi"
#	
#	1. Set machine with enviroment variable
#
#	MACHINE = set in ~/.bashrc or ~/.bashrc.ext
#
#

ifeq ($(MACHINE),dlghp)
	HOME = /home/dg6/code
	USEGPU=yes
	ifeq ($(USEGPU),yes)
		CU = /home/dg6/code/cuda/bin/nvcc
		CUFLAGS = -g -O3 -arch=sm_11 --ptxas-options=-v --compiler-bindir /home/dg6/code/gcc43/usr/bin -I /usr/include/c++/4.4.0/x86_64-redhat-linux/ -I /usr/include/c++/4.4.0/
		CUDA_LIB = -L/home/dg6/code/cuda/lib/ -lcudart -fopenmp
		CPP = g++
		CPPFLAGS = -g -fbacktrace
		CUDA_DIR = src/cuda
	endif
	CPP_DIRECTIVES = -DUSE_DISLIN=1 -DLGHP
endif

ifeq ($(MACHINE),lens)
	USEGPU=yes
	ifeq ($(USEGPU),yes)
		CU = nvcc
		CUFLAGS = -g -O3 -arch=sm_13 --ptxas-options=-v 
		CUDA_LIB = -L /sw/analysis-x64/cuda/2.2/sl5.0_binary/lib -lcudart -mp 
		CPP = g++
		CPPFLAGS = -g 
		CUDA_DIR = src/cuda
	endif
	CPP_DIRECTIVES = -DUSE_DISLIN=0 -DLGHP
endif

ifeq ($(MACHINE),$(filter $(MACHINE),franklin jaguar))
	CRAYXT4=1
endif

EXEC = ${PWD}/xaorsa2d.$(MACHINE)

SRC_DIR = ${PWD}/src

OBJ_DIR = ${PWD}/obj

MOD_DIR = ${PWD}/mod

FFT_DIR = $(SRC_DIR)/fftpack

CQL_DIR = $(SRC_DIR)/cql3d

INCLUDE_DIR = $(SRC_DIR)/cql3d

DLG_DIR = $(SRC_DIR)/dlg

OBJ_FILES = \
 $(OBJ_DIR)/cauchy_mod.o \
 $(OBJ_DIR)/size_mod.o \
 $(OBJ_DIR)/aorsa2din_mod.o \
 $(OBJ_DIR)/swim_global_data_mod.o \
 $(OBJ_DIR)/precision_mod.o \
 $(OBJ_DIR)/profile_mod.o\
 $(OBJ_DIR)/qlsum.o \
 $(OBJ_DIR)/ql_myra.o \
 $(OBJ_DIR)/mets2aorsa.o \
 $(OBJ_DIR)/cauchy_ppart.o \
 $(OBJ_DIR)/vlog.o \
 $(OBJ_DIR)/aorsaSubs.o \
 $(OBJ_DIR)/sigma.o \
 $(OBJ_DIR)/zfunction.o \
 $(OBJ_DIR)/ztable.o \
 $(OBJ_DIR)/current.o \
 $(OBJ_DIR)/mets2aorsa_myra.o \
 $(OBJ_DIR)/slowDown.o \
 $(OBJ_DIR)/fourier.o \
 $(OBJ_DIR)/assert.o \
 $(OBJ_DIR)/setupblacs.o \
 $(OBJ_DIR)/bessel.o \
 $(OBJ_DIR)/check.o \
 $(OBJ_DIR)/rf2x_setup2.o \
 $(OBJ_DIR)/profile_setup.o \
 $(OBJ_DIR)/eqdsk_setup.o \
 $(OBJ_DIR)/orbit.o \
 $(OBJ_DIR)/eqdsk_plot.o \
 $(OBJ_DIR)/fieldws.o \
 $(OBJ_DIR)/ranlib.o \
 $(OBJ_DIR)/dshell.o \
 $(OBJ_DIR)/aorsa2dMain.o 

ifeq ($(USEGPU),yes)
OBJ_CUDA = \
		   $(OBJ_DIR)/qlsum_gpu_kernels.o \
		   $(OBJ_DIR)/qlsum_gpu_host.o 
endif

OBJ_FFT = \
 $(OBJ_DIR)/cfftb1.o \
 $(OBJ_DIR)/cfftf1.o \
 $(OBJ_DIR)/cffti1.o \
 $(OBJ_DIR)/passb.o \
 $(OBJ_DIR)/passb2.o \
 $(OBJ_DIR)/passb3.o \
 $(OBJ_DIR)/passb4.o \
 $(OBJ_DIR)/passb5.o \
 $(OBJ_DIR)/passf.o \
 $(OBJ_DIR)/passf2.o \
 $(OBJ_DIR)/passf3.o \
 $(OBJ_DIR)/passf4.o \
 $(OBJ_DIR)/passf5.o \
 $(OBJ_DIR)/zfftb.o \
 $(OBJ_DIR)/zfftf.o \
 $(OBJ_DIR)/zffti.o
 
OBJ_CQL3D_SETUP = \
 $(OBJ_DIR)/global_data_m.o \
 $(OBJ_DIR)/basis_functions_m.o \
 $(OBJ_DIR)/f_expanded_m.o \
 $(OBJ_DIR)/CQL_kinds_m.o \
 $(OBJ_DIR)/vector_write_m.o \
 $(OBJ_DIR)/read_cql3d.o \
 $(OBJ_DIR)/fitpack.o \
 $(OBJ_DIR)/cubic_B_splines_v.o \
 $(OBJ_DIR)/cql3d_setup.o

OBJ_DLG = \
 $(OBJ_DIR)/read_particle_f.o \
 $(OBJ_DIR)/gc_integrate.o \
 $(OBJ_DIR)/interp.o \
 $(OBJ_DIR)/eqdsk_dlg.o \
 $(OBJ_DIR)/gc_terms.o \
 $(OBJ_DIR)/constants.o \
 $(OBJ_DIR)/dlg.o \
 ${OBJ_DIR}/write_pql.o \
 ${OBJ_DIR}/write_pf.o


ifeq ($(CRAYXT4),1)
	OPT_FLAGS = #-fast -tp barcelona-64 -traceback #-Minfo -mp#-Mlist #-Mconcur=nonuma
	PRECISION_FLAGS = -Mr8 -r8
	COMMON_OPTION  = -Msave -I $(INCLUDE_DIR) $(OPT_FLAGS) -module $(MOD_DIR) -Kieee
	COMMON_OPTION2  = -Mnosave -I $(INCLUDE_DIR) $(OPT_FLAGS) -module $(MOD_DIR) -Kieee
	F77 = ftn $(COMMON_OPTION) $(PRECISION_FLAGS)
	F77_NOSAVE = ftn $(COMMON_OPTION2) $(PRECISION_FLAGS)
	F77_r4= ftn $(COMMON_OPTION2)
	F77_FFT = ftn $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90 = ftn -Mfree $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90_lst = ftn -Mfree $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90_NOSAVE =   ftn -Mfree $(COMMON_OPTION2) $(PRECISION_FLAGS)
	F90_2 = ftn -Mfixed $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90_2_NOSAVE = ftn -Mfixed $(COMMON_OPTION2) $(PRECISION_FLAGS)
	F90_CQL = ftn -Mfree $(COMMON_OPTION2) $(PRECISION_FLAGS)
	COMPILE_DLG	= ftn -c -I $(INCLUDE_DIR) $(PRECISION_FLAGS) $(OPT_FLAGS) -module $(MOD_DIR)
	COMPILE_DLG_R4	= ftn -c -I $(INCLUDE_DIR) $(OPT_FLAGS) -module $(MOD_DIR)
	OPTIMIZATION = 
	F77FLAGS= -c $(OPTIMIZATION)  -I $(INCLUDE_DIR) 
	BOUNDS = 
	WARN =

	INC_DIR = 	
	BOUNDS = # turning -Mbounds on causes PBLAS to crash
	WARN = -Minform=warn
	DEBUG = -g -traceback
	FFLAGS = ${WARN} ${DEBUG} 
	F90FLAGS = ${WARN} ${DEBUG}
	FFTFLAGS = -Msave

endif

ifeq ($(MACHINE),dlghp)

	MPI_INCLUDE_DIR = $(HOME)/openmpi/gnu_64/include

	MPIF90 = $(HOME)/openmpi/gnu_64/bin/mpif90
	MPIF77 = $(HOME)/openmpi/gnu_64/bin/mpif77

	COMMON_OPTION	= -fdefault-real-8 -fno-automatic -I $(MPI_INCLUDE_DIR) -J$(MOD_DIR)
	COMMON_OPTION2	= -fdefault-real-8 -fautomatic -I $(MPI_INCLUDE_DIR) -J$(MOD_DIR)
	COMMON_OPTION3	= -fautomatic -I -I $(MPI_INCLUDE_DIR) -J$(MOD_DIR)
	#-finit-local-zero 
	F77	= $(MPIF77)  $(COMMON_OPTION)
	F77_NOSAVE	= $(MPIF77) $(COMMON_OPTION2)
	F77_r4	= $(MPIF77) $(COMMON_OPTION3)
	F77_FFT	= $(MPIF77) $(COMMON_OPTION)
	F90	= $(MPIF90) -ffree-form $(COMMON_OPTION)
	F90_NOSAVE	= $(MPIF90) -ffree-form $(COMMON_OPTION2)
	F90_2 = $(MPIF90) -ffixed-form $(COMMON_OPTION)
	F90_2_NOSAVE = $(MPIF90) -ffixed-form $(COMMON_OPTION2)
	F90_CQL	= $(MPIF90) -ffree-form $(COMMON_OPTION2)
	COMPILE_DLG	= $(MPIF90) -c -I $(INCLUDE_DIR) -I $(MPI_INCLUDE_DIR) -I $(HOME)/netcdf/netcdf_gnu64/include/ -fdefault-real-8 -J$(MOD_DIR)
	COMPILE_DLG_R4	= $(MPIF77) -c -I $(INCLUDE_DIR) -I $(MPI_INCLUDE_DIR) -I $(HOME)/netcdf/netcdf_gnu64/include/ -J$(MOD_DIR)
	OPTIMIZATION = #-O3 -march=core2 -funroll-loops -Wuninitialized
	F77FLAGS= -c $(OPTIMIZATION)  -I $(INCLUDE_DIR) -I$(HOME)/netcdf/netcdf_gnu64/include/
	F90_lst = $(MPIF90) -ffree-form $(COMMON_OPTION)

	INC_DIR = 	
	BOUNDS = -fbounds-check
	WARN = -Wall
	DEBUG = -g -fbacktrace
	FFLAGS = ${WARN} ${DEBUG} 
	F90FLAGS = ${WARN} ${DEBUG}

endif
	
ifeq ($(MACHINE),lens)
	MPIF77 = mpif90 -module $(MOD_DIR) -traceback
	MPIF90 = mpif90 -module $(MOD_DIR) -traceback
	OPT_FLAGS = #-fast #-tp barcelona-64 -traceback #-Minfo -mp#-Mlist #-Mconcur=nonuma
	PRECISION_FLAGS = -Mr8 -r8
	COMMON_OPTION  = -Msave -I $(INCLUDE_DIR) $(OPT_FLAGS)  -Kieee
	COMMON_OPTION2  = -Mnosave -I $(INCLUDE_DIR) $(OPT_FLAGS) -Kieee
	F77 = ${MPIF77} $(COMMON_OPTION) $(PRECISION_FLAGS)
	F77_NOSAVE = ${MPIF77} $(COMMON_OPTION2) $(PRECISION_FLAGS)
	F77_r4= ${MPIF77} $(COMMON_OPTION2)
	F77_FFT = ${MPIF77} $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90 = ${MPIF90} -Mfree $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90_lst = ${MPIF90} -Mfree $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90_NOSAVE =   ${MPIF90} -Mfree $(COMMON_OPTION2) $(PRECISION_FLAGS)
	F90_2 = ${MPIF90} -Mfixed $(COMMON_OPTION) $(PRECISION_FLAGS)
	F90_2_NOSAVE = ${MPIF90} -Mfixed $(COMMON_OPTION2) $(PRECISION_FLAGS)
	F90_CQL = ${MPIF90} -Mfree $(COMMON_OPTION2) $(PRECISION_FLAGS)
	COMPILE_DLG	= ${MPIF90} -c -I $(INCLUDE_DIR) $(PRECISION_FLAGS) $(OPT_FLAGS) -module $(MOD_DIR)
	COMPILE_DLG_R4	= ${MPIF90} -c -I $(INCLUDE_DIR) $(OPT_FLAGS) -module $(MOD_DIR)
	OPTIMIZATION = 
	F77FLAGS= -c $(OPTIMIZATION)  -I $(INCLUDE_DIR) 
	BOUNDS = 
	WARN =
endif


INLINE=


OPTIMIZATION4 = $(OPTIMIZATION)
OPTIMIZATIONg =  

F77FLAGSg=  -c  -R abcps \
	$(OPTIMIZATIONg) -I $(INCLUDE_DIR)
F77FLAGS4= -c $(OPTIMIZATION4) -I $(INCLUDE_DIR)

COMPILE = $(F77) $(F77FLAGS)
COMPILE_NOSAVE = $(F77_NOSAVE) $(F77FLAGS)
COMPILE_r4 = $(F77_r4) $(F77FLAGS)

COMPILE90 = $(F90) $(F77FLAGS) $(INLINE)
COMPILE90_lst = $(F90_lst) $(F77FLAGS) $(INLINE)

COMPILE90_NOSAVE = $(F90_NOSAVE) $(F77FLAGS) $(INLINE)
COMPILE90_2 = $(F90_2) $(F77FLAGS) $(INLINE)
COMPILE90_2_NOSAVE = $(F90_2_NOSAVE) $(F77FLAGS) $(INLINE)
COMPILEg = $(F77) $(F77FLAGSg)
COMPILE4 = $(F77) $(F77FLAGS4)

COMPILE_FFT  = $(F77_FFT)  $(F77FLAGS)

COMPILE_CQL  = $(F90_CQL)  $(F77FLAGS)

#HPL_LIB = /ccs/home/efdazedo/public/jaguar/libzhpl.a


ifeq ($(MACHINE),dlghp)

	DISLIN_DIR = /home/dg6/code/dislin/dislin64
	DISLIN = -I ${DISLIN_DIR}/gf -L ${DISLIN_DIR} -ldislin
	BLACS = \
		$(HOME)/blacs/blacs_gnu64/LIB/blacs_MPI-LINUX-0.a \
		$(HOME)/blacs/blacs_gnu64/LIB/blacsCinit_MPI-LINUX-0.a \
		$(HOME)/blacs/blacs_gnu64/LIB/blacsF77init_MPI-LINUX-0.a \
		$(HOME)/blacs/blacs_gnu64/LIB/blacs_MPI-LINUX-0.a
	BLAS = $(HOME)/goto_blas/GotoBLAS/libgoto.a -lpthread
	HPL = $(HOME)/hpl/hpl-2.0/lib/linux/libhpl.a
	FFT_LIB = $(HOME)/fftpack/dfftpack_gnu64/libdfftpack.a
	SCALAPACK = $(HOME)/scalapack/scalapack_gnu64/libscalapack.a
	X11 = -L/usr/lib -lX11
	LOADFLAGS = 
	PNETCDF_DIR = $(HOME)/pNetCdf/pnetcdf_gnu64
	PNETCDF = -I ${PNETCDF_DIR}/include ${PNETCDF_DIR}/lib/libpnetcdf.a
	NETCDF_DIR = $(HOME)/netcdf/netcdf_gnu64
	NETCDF = -I ${NETCDF_DIR}/include
	LIBS = -L$(HOME)/netcdf/netcdf_gnu64/lib -lnetcdf $(HOME)/pgplot/pgplot_gnu64/libpgplot.a $(SCALAPACK) $(HPL) $(BLACS) $(BLAS) $(X11) $(FFT_LIB) ${PNETCDF}
	INC_DIR = -I ${NETCDF_DIR}/include -I ${DISLIN_DIR}/gf -I ${PNETCDF_DIR}/include
endif

ifeq ($(MACHINE),franklin)
	HPL_LIB = /u0/j/jaegeref/AORSA2D/libzhpl.a
	FFT_LIB = /u0/j/jaegeref/AORSA2D/libdfftpack_pgf90.a
	LOADFLAGS = -lsci
	LIBS = \
		$(HPL_LIB) \
		$(FFT_LIB) \
	  	${NETCDF} \
		${PNETCDF} \
		-L/u0/j/jaegeref/pgplot_xt3 -lpgplot #-L/ccs/home/jaegeref/pgplot_xt3 -lpgplot 
endif

ifeq ($(MACHINE),jaguar)
	HPL_LIB = /ccs/home/jaegeref/AORSA2D/libzhpl.a
	FFT_LIB = /ccs/home/jaegeref/AORSA2D/libdfftpack_pgf90.a
	LOADFLAGS =

	PNETCDF = ${PNETCDF_LIB}
	NETCDF = ${NETCDF_FLIB}
	LIBS = \
	  $(HPL_LIB) \
	  $(FFT_LIB) \
	  ${PNETCDF} \
	  ${NETCDF_FLIB} $(PGPLOT_LIB) 
	INC_DIR = -I ${NETCDF_DIR}/include -I ${PNETCDF_DIR}/include

endif

ifeq ($(MACHINE),lens)
	#SCI_LIB = -L bbales/lib/ -lscalapack -lblacs -lblacsC -lblacsF77 -lblacs 
	FFT_LIB = /ccs/home/jaegeref/AORSA2D/libdfftpack_pgf90.a
	LOADFLAGS =

	PNETCDF = 
	NETCDF = ${NETCDF_FLIB}
	LIBS = \
	  	$(SCI_LIB) \
	  	$(FFT_LIB) \
	  	${PNETCDF_LIB} \
	  	${NETCDF_FLIB} \
		$(PGPLOT_LIB) \
		${SCALAPACK_LIB} \
		${BLACS_LIB} \
		${BLASGOTO_LIB} 
	INC_DIR = -I ${NETCDF_DIR}/include -I ${PNETCDF_DIR}/include
endif


PERFORMANCE = 
#-L/ccs/home/larkin/xt3/fpmpi/lib -lfpmpi_papi

LOAD = $(F90) $(OPTIMIZATION) $(LOADFLAGS) 

# Compile the program

.PHONY: depend clean

ifeq ($(MACHINE),dlghp)
$(EXEC): $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) ${OBJ_CUDA}
	$(LOAD) -o $(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) $(LIBS) $(DISLIN) ${OBJ_CUDA} ${CUDA_LIB} 
else
$(EXEC): load_modules $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) ${OBJ_CUDA}
	$(LOAD) -o $(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) $(LIBS) ${OBJ_CUDA} ${CUDA_LIB}
endif

load_modules:
#	${module} load ${NETCDF_MODULE} 
#	${module} load ${PNETCDF_MODULE}
#	${module} load ${PGPLOT_MODULE} 

# Dependencies

ifeq ($(USEGPU),yes)

$(OBJ_DIR)/%.o : $(CUDA_DIR)/%.cu
	$(CU) -c $? $(CUFLAGS) -o $@

$(OBJ_DIR)/%.o : $(CUDA_DIR)/%.cpp
	$(CPP) -c $? $(CPPFLAGS) -o $@

endif

$(OBJ_DIR)/mets2aorsa.o:     $(SRC_DIR)/mets2aorsa.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/mets2aorsa.o \
                             $(SRC_DIR)/mets2aorsa.f ${BOUNDS}
				
$(OBJ_DIR)/mets2aorsa_myra.o: $(SRC_DIR)/mets2aorsa_myra.f 
	                      $(COMPILE90) -o $(OBJ_DIR)/mets2aorsa_myra.o \
                              $(SRC_DIR)/mets2aorsa_myra.f ${BOUNDS}

ifeq ($(USEGPU),yes)
$(OBJ_DIR)/qlsum.o:  $(CUDA_DIR)/qlsum_gpu.f90
	$(COMPILE_DLG) -o $(OBJ_DIR)/qlsum.o $(CUDA_DIR)/qlsum_gpu.f90 $(NETCDF) ${BOUNDS} ${WARN}
else
$(OBJ_DIR)/qlsum.o:  $(SRC_DIR)/qlsum.f90 
	$(COMPILE_DLG) -o $(OBJ_DIR)/qlsum.o $(SRC_DIR)/qlsum.f90 $(NETCDF) ${BOUNDS} ${WARN}
endif			

$(OBJ_DIR)/cauchy_ppart.o:   $(SRC_DIR)/cauchy_ppart.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/cauchy_ppart.o \
                             $(SRC_DIR)/cauchy_ppart.f ${BOUNDS}	

$(OBJ_DIR)/vlog.o:           $(SRC_DIR)/vlog.f 
		             $(COMPILE) -o $(OBJ_DIR)/vlog.o \
			     $(SRC_DIR)/vlog.f ${BOUNDS}
			     
$(OBJ_DIR)/aorsa2dMain.o:	$(SRC_DIR)/aorsa2dMain.F90 $(OBJ_DIR)/interp.o ${OBJ_DIR}/write_pf.o 
	$(COMPILE_DLG) ${CPP_DIRECTIVES} -o $(OBJ_DIR)/aorsa2dMain.o $(SRC_DIR)/aorsa2dMain.F90 $(NETCDF) ${BOUNDS}			     

$(OBJ_DIR)/dshell.o:         $(SRC_DIR)/dshell.f 
	                     $(COMPILE90_2) -o $(OBJ_DIR)/dshell.o \
                             $(SRC_DIR)/dshell.f ${BOUNDS}

$(OBJ_DIR)/aorsaSubs.o:      $(SRC_DIR)/aorsaSubs.f 
	                     $(COMPILE90_2) -o $(OBJ_DIR)/aorsaSubs.o \
                             $(SRC_DIR)/aorsaSubs.f ${BOUNDS}

$(OBJ_DIR)/sigma.o:          $(SRC_DIR)/sigma.f  $(OBJ_DIR)/read_particle_f.o ${OBJ_DIR}/write_pf.o
	                     $(COMPILE90_2) -o $(OBJ_DIR)/sigma.o \
                             $(SRC_DIR)/sigma.f ${NETCDF} ${BOUNDS}
                                
$(OBJ_DIR)/zfunction.o:      $(SRC_DIR)/zfunction.f 
	                     $(COMPILE90_2) -o $(OBJ_DIR)/zfunction.o \
                             $(SRC_DIR)/zfunction.f ${BOUNDS}
                                
$(OBJ_DIR)/ztable.o:         $(SRC_DIR)/ztable.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/ztable.o \
                             $(SRC_DIR)/ztable.f ${BOUNDS}							                                 
                                 
$(OBJ_DIR)/bessel.o:         $(SRC_DIR)/bessel.f 
	                     $(COMPILE) -o $(OBJ_DIR)/bessel.o \
                             $(SRC_DIR)/bessel.f ${BOUNDS}                                
                                                              
$(OBJ_DIR)/current.o:        $(SRC_DIR)/current.f 
	                     $(COMPILE90_2) -o $(OBJ_DIR)/current.o \
                             $(SRC_DIR)/current.f ${BOUNDS}

ifeq ($(USEGPU),yes)
$(OBJ_DIR)/ql_myra.o:	$(CUDA_DIR)/ql_myra_gpu.f $(OBJ_DIR)/read_particle_f.o $(OBJ_DIR)/gc_integrate.o ${OBJ_DIR}/write_pql.o 
	$(COMPILE_NOSAVE) ${CPP_DIRECTIVES} -o $(OBJ_DIR)/ql_myra.o $(CUDA_DIR)/ql_myra_gpu.f $(NETCDF) ${BOUNDS} ${PNETCDF}
else
$(OBJ_DIR)/ql_myra.o:	$(SRC_DIR)/ql_myra.F $(OBJ_DIR)/read_particle_f.o $(OBJ_DIR)/gc_integrate.o ${OBJ_DIR}/write_pql.o 
	$(COMPILE_NOSAVE) ${CPP_DIRECTIVES} -o $(OBJ_DIR)/ql_myra.o $(SRC_DIR)/ql_myra.F $(NETCDF) ${BOUNDS} ${PNETCDF}
endif

$(OBJ_DIR)/slowDown.o:       $(SRC_DIR)/slowDown.f 
	                     $(COMPILE) -o $(OBJ_DIR)/slowDown.o \
                             $(SRC_DIR)/slowDown.f ${BOUNDS}                                                                           

$(OBJ_DIR)/fourier.o:        $(SRC_DIR)/fourier.f 
	                     $(COMPILE90_2) -o $(OBJ_DIR)/fourier.o \
                             $(SRC_DIR)/fourier.f ${BOUNDS}
                                
$(OBJ_DIR)/assert.o:         $(SRC_DIR)/assert.f 
	                     $(COMPILE) -o $(OBJ_DIR)/assert.o \
                             $(SRC_DIR)/assert.f ${BOUNDS}                            
                                
$(OBJ_DIR)/setupblacs.o:     $(SRC_DIR)/setupblacs.f 
	                     $(COMPILE) -o $(OBJ_DIR)/setupblacs.o \
                             $(SRC_DIR)/setupblacs.f ${BOUNDS}                                                                                                

$(OBJ_DIR)/check.o:          $(SRC_DIR)/check.f 
	                     $(COMPILE) -o $(OBJ_DIR)/check.o \
                             $(SRC_DIR)/check.f ${BOUNDS}
			     
$(OBJ_DIR)/cauchy_mod.o:    $(SRC_DIR)/cauchy_mod.f 
	                     $(COMPILE) -o $(OBJ_DIR)/cauchy_mod.o \
                             $(SRC_DIR)/cauchy_mod.f ${BOUNDS}			     
				
$(OBJ_DIR)/precision_mod.o:  $(SRC_DIR)/precision_mod.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/precision_mod.o \
                             $(SRC_DIR)/precision_mod.f ${BOUNDS}
ifeq ($(MACHINE),dlghp)			     
$(OBJ_DIR)/size_mod.o:       $(SRC_DIR)/size_mod.F 
	                     $(COMPILE90) -DDLGHP=1 -o $(OBJ_DIR)/size_mod.o \
                             $(SRC_DIR)/size_mod.F ${BOUNDS}			     			     
else
$(OBJ_DIR)/size_mod.o:       $(SRC_DIR)/size_mod.F 
	                     $(COMPILE90) -DDLGHP=0 -o $(OBJ_DIR)/size_mod.o \
                             $(SRC_DIR)/size_mod.F ${BOUNDS}			     			     
endif
		
$(OBJ_DIR)/aorsa2din_mod.o:  $(SRC_DIR)/aorsa2din_mod.f $(INC_FILES)
	                     $(COMPILE90) -o $(OBJ_DIR)/aorsa2din_mod.o \
                             $(SRC_DIR)/aorsa2din_mod.f ${BOUNDS}					 

$(OBJ_DIR)/profile_mod.o:    $(SRC_DIR)/profile_mod.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/profile_mod.o \
                             $(SRC_DIR)/profile_mod.f ${BOUNDS}
			     
$(OBJ_DIR)/swim_global_data_mod.o: $(SRC_DIR)/swim_global_data_mod.f 
	                           $(COMPILE90) -o $(OBJ_DIR)/swim_global_data_mod.o \
                                   $(SRC_DIR)/swim_global_data_mod.f ${BOUNDS}
		
$(OBJ_DIR)/rf2x_setup2.o:    $(SRC_DIR)/rf2x_setup2.f 
	                     $(COMPILE90_2_NOSAVE) -o $(OBJ_DIR)/rf2x_setup2.o \
                             $(SRC_DIR)/rf2x_setup2.f ${BOUNDS}
			     
$(OBJ_DIR)/profile_setup.o:  $(SRC_DIR)/profile_setup.f 
	                     $(COMPILE90_2_NOSAVE) -o $(OBJ_DIR)/profile_setup.o \
                             $(SRC_DIR)/profile_setup.f 
			     
$(OBJ_DIR)/eqdsk_setup.o:    $(SRC_DIR)/eqdsk_setup.f90 
	                     $(COMPILE_DLG) -o $(OBJ_DIR)/eqdsk_setup.o \
                             $(SRC_DIR)/eqdsk_setup.f90 ${BOUNDS}
				
$(OBJ_DIR)/orbit.o:          $(SRC_DIR)/orbit.f 
	                     $(COMPILE90_2_NOSAVE) -o $(OBJ_DIR)/orbit.o \
                             $(SRC_DIR)/orbit.f			     
			     			     			     			     
$(OBJ_DIR)/eqdsk_plot.o:     $(SRC_DIR)/eqdsk_plot.f90 ${OBJ_DIR}/fieldws.o
	                     $(COMPILE_DLG_R4) -o $(OBJ_DIR)/eqdsk_plot.o \
                             $(SRC_DIR)/eqdsk_plot.f90 ${BOUNDS}				     
			     			     
$(OBJ_DIR)/fieldws.o:        $(SRC_DIR)/fieldws.f90
	                     $(COMPILE_DLG_R4) -o $(OBJ_DIR)/fieldws.o \
                             $(SRC_DIR)/fieldws.f90 $(NETCDF) ${BOUNDS}

			    		    		     		   			     			     				
# FFT files
# ---------

${OBJ_DIR}/%.o: ${FFT_DIR}/%.f
	${F77} -c ${FFLAGS} $< -o $@ ${FFTFLAGS} 
				
# SRC files
# ---------

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f
	${F77} -c ${FFLAGS} $< -o $@ ${BOUNDS} ${NETCDF}

# CQL files
# ---------

${OBJ_DIR}/%.o: ${CQL_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${NETCDF}

${OBJ_DIR}/%.o: ${CQL_DIR}/%.f
	${F77} -c ${FFLAGS} $< -o $@ ${BOUNDS} ${NETCDF}

# DLG files
# ---------
		
${OBJ_DIR}/%.o: ${DLG_DIR}/%.F
	${F77} -c ${FFLAGS} $< -o $@ ${INC_DIR} ${CPP_DIRECTIVES}

${OBJ_DIR}/%.o: ${DLG_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${INC_DIR} ${CPP_DIRECTIVES}

${OBJ_DIR}/%.o: ${DLG_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${INC_DIR}

${OBJ_DIR}/%.o: ${DLG_DIR}/%.f
	${F77} -c ${FFLAGS} $< -o $@

clean:
	rm $(EXEC) $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

MAKEFILE_INC=Makefile.deps
SRCDIRS=${DLG_DIR} ${CQL_DIR} ${SRC_DIR}
FSRCS0:=$(foreach DIR, . $(SRCDIRS),$(wildcard $(DIR)/*.f90 $(DIR)/*.F ${DIR}/*.f ${DIR}/*.F90))
FSRCS:=$(sort $(notdir $(FSRCS0)))
Includes=#-I${PNETCDF_DIR}/include

F_makedepend=./sfmakedepend --file - $(addprefix --srcdir ,$(SRCDIRS)) $(subst -I,-I ,$(Includes)) --objdir ${OBJ_DIR} --moddir ${OBJ_DIR} --depend=obj
depend $(MAKEFILE_INC): 
	$(F_makedepend) $(FSRCS) > $(MAKEFILE_INC)

include $(MAKEFILE_INC)

#depend: 
#	./sfmakedepend --depend=obj --objdir ${OBJ_DIR} --moddir ${OBJ_DIR} ${DLG_DIR}/*.f90 



