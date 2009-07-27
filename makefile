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
#	2. Also set CODEDIR below
#

ifeq ($(MACHINE),dlghp)
	HOME = /home/dg6/code
	DISLIN = -I /home/dg6/code/dislin/dislin64/gf -L /home/dg6/code/dislin/dislin64 -ldislin
endif
ifeq ($(MACHINE),franklin)
	CRAYXT4=1
	NETCDF_MODULE = netcdf
	PNETCDF_MODULE = p-netcdf
	PGPLOT_MODULE = 
endif
ifeq ($(MACHINE),jaguar)
	CRAYXT4=1
	NETCDF_MODULE = netcdf
	PNETCDF_MODULE = p-netcdf
	PGPLOT_MODULE = pgplot
	module = /opt/modules/default/bin/modulecmd bash
endif

CODEDIR = AORSA2D/AORSA-17-DLG

EXEC = xaorsa2d.$(MACHINE)

SRC_DIR = $(HOME)/$(CODEDIR)/src

OBJ_DIR = $(HOME)/$(CODEDIR)/obj/$(MACHINE)

MOD_DIR = $(HOME)/$(CODEDIR)/mod/$(MACHINE)

FFT_DIR = $(SRC_DIR)/fftpack

CQL3D_SETUP_DIR = $(SRC_DIR)/cql3d

INCLUDE_DIR = $(SRC_DIR)/CQL3D_SETUP

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
endif

ifeq ($(MACHINE),dlghp)
	MPI_INCLUDE_DIR = $(HOME)/openmpi/gnu_64/include
	MPIF90 = $(HOME)/openmpi/gnu_64/bin/mpif90 -g -fbacktrace
	MPIF77 = $(HOME)/openmpi/gnu_64/bin/mpif77 -g -fbacktrace
	COMMON_OPTION	= -fdefault-real-8 -fno-automatic  -I $(INCLUDE_DIR) -I $(MPI_INCLUDE_DIR) -J$(MOD_DIR)
	COMMON_OPTION2	= -fdefault-real-8 -fautomatic -I $(INCLUDE_DIR) -I $(MPI_INCLUDE_DIR) -J$(MOD_DIR)
	COMMON_OPTION3	= -fautomatic -I $(INCLUDE_DIR) -I $(MPI_INCLUDE_DIR) -J$(MOD_DIR)
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
	BOUNDS = -fbounds-check
	WARN = -Wall
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
	PNETCDF = -I $(HOME)/pNetCdf/pnetcdf_gnu64/include $(HOME)/pNetCdf/pnetcdf_gnu64/lib/libpnetcdf.a
	NETCDF = -I $(HOME)/netcdf/netcdf_gnu64/lib
	LIBS = -L$(HOME)/netcdf/netcdf_gnu64/lib -lnetcdf $(HOME)/pgplot/pgplot_gnu64/libpgplot.a $(SCALAPACK) $(HPL) $(BLACS) $(BLAS) $(X11) $(FFT_LIB) ${PNETCDF}
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
endif

PERFORMANCE = 
#-L/ccs/home/larkin/xt3/fpmpi/lib -lfpmpi_papi

LOAD = $(F90) $(OPTIMIZATION) $(LOADFLAGS) 

# Compile the program

ifeq ($(MACHINE),dlghp)
$(EXEC): $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG)
	$(LOAD) -o ./$(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) $(LIBS) $(DISLIN)
else
$(EXEC): load_modules $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG)
	$(LOAD) -o ./$(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(OBJ_DLG) $(LIBS) 
endif

load_modules:
#	${module} load ${NETCDF_MODULE} 
#	${module} load ${PNETCDF_MODULE}
#	${module} load ${PGPLOT_MODULE} 

# Dependencies


$(OBJ_DIR)/mets2aorsa.o:     $(SRC_DIR)/mets2aorsa.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/mets2aorsa.o \
                             $(SRC_DIR)/mets2aorsa.f ${BOUNDS}
				
$(OBJ_DIR)/mets2aorsa_myra.o: $(SRC_DIR)/mets2aorsa_myra.f 
	                      $(COMPILE90) -o $(OBJ_DIR)/mets2aorsa_myra.o \
                              $(SRC_DIR)/mets2aorsa_myra.f ${BOUNDS}
				 
$(OBJ_DIR)/qlsum.o:  $(SRC_DIR)/qlsum.f90
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/qlsum.o \
                              $(SRC_DIR)/qlsum.f90 $(NETCDF) ${BOUNDS} ${WARN}
				
$(OBJ_DIR)/cauchy_ppart.o:   $(SRC_DIR)/cauchy_ppart.f 
	                     $(COMPILE90) -o $(OBJ_DIR)/cauchy_ppart.o \
                             $(SRC_DIR)/cauchy_ppart.f ${BOUNDS}	

$(OBJ_DIR)/vlog.o:           $(SRC_DIR)/vlog.f 
		             $(COMPILE) -o $(OBJ_DIR)/vlog.o \
			     $(SRC_DIR)/vlog.f ${BOUNDS}
			     
ifeq ($(MACHINE),dlghp)
$(OBJ_DIR)/aorsa2dMain.o:    $(SRC_DIR)/aorsa2dMain.F90 $(OBJ_DIR)/interp.o ${OBJ_DIR}/write_pf.o
	                     $(COMPILE_DLG) -o $(OBJ_DIR)/aorsa2dMain.o \
                             $(SRC_DIR)/aorsa2dMain.F90 $(NETCDF) ${BOUNDS}			     
else
$(OBJ_DIR)/aorsa2dMain.o:    $(SRC_DIR)/aorsa2dMain.F90 $(OBJ_DIR)/interp.o ${OBJ_DIR}/write_pf.o
	                     $(COMPILE_DLG) -DUSE_HPL -o $(OBJ_DIR)/aorsa2dMain.o \
                             $(SRC_DIR)/aorsa2dMain.F90 $(NETCDF) ${BOUNDS}		     
endif

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

ifeq (${MACHINE},dlghp)
$(OBJ_DIR)/ql_myra.o:        $(SRC_DIR)/ql_myra.F $(OBJ_DIR)/read_particle_f.o $(OBJ_DIR)/gc_integrate.o ${OBJ_DIR}/write_pql.o
	                     $(COMPILE_NOSAVE) -DDLGHP -o $(OBJ_DIR)/ql_myra.o \
                             $(SRC_DIR)/ql_myra.F $(NETCDF) ${BOUNDS} ${PNETCDF}
else
$(OBJ_DIR)/ql_myra.o:        $(SRC_DIR)/ql_myra.F $(OBJ_DIR)/read_particle_f.o $(OBJ_DIR)/gc_integrate.o ${OBJ_DIR}/write_pql.o
	                     $(COMPILE_NOSAVE) -DDLGHP -o $(OBJ_DIR)/ql_myra.o \
                             $(SRC_DIR)/ql_myra.F $(NETCDF) ${BOUNDS} ${PNETCDF}
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
			     			     			     			     
$(OBJ_DIR)/eqdsk_plot.o:     $(SRC_DIR)/eqdsk_plot.f 
	                     $(COMPILE_r4) -o $(OBJ_DIR)/eqdsk_plot.o \
                             $(SRC_DIR)/eqdsk_plot.f ${BOUNDS}				     
			     			     
$(OBJ_DIR)/fieldws.o:        $(SRC_DIR)/fieldws.f90
	                     $(COMPILE_DLG_R4) -o $(OBJ_DIR)/fieldws.o \
                             $(SRC_DIR)/fieldws.f90 $(NETCDF) ${BOUNDS}

		     			     
			     
			    		    		     		   			     			     				
### FFTPACK files:

$(OBJ_DIR)/cfftb1.o :      $(FFT_DIR)/cfftb1.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/cfftb1.o \
			   $(FFT_DIR)/cfftb1.f ${BOUNDS}

$(OBJ_DIR)/cfftf1.o :      $(FFT_DIR)/cfftf1.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/cfftf1.o \
			   $(FFT_DIR)/cfftf1.f ${BOUNDS}

$(OBJ_DIR)/cffti1.o:       $(FFT_DIR)/cffti1.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/cffti1.o \
			   $(FFT_DIR)/cffti1.f ${BOUNDS}

$(OBJ_DIR)/passb.o :       $(FFT_DIR)/passb.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passb.o \
			   $(FFT_DIR)/passb.f ${BOUNDS}

$(OBJ_DIR)/passb2.o :      $(FFT_DIR)/passb2.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passb2.o \
			   $(FFT_DIR)/passb2.f 

$(OBJ_DIR)/passb3.o :      $(FFT_DIR)/passb3.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passb3.o \
			   $(FFT_DIR)/passb3.f ${BOUNDS}

$(OBJ_DIR)/passb4.o :      $(FFT_DIR)/passb4.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passb4.o \
			   $(FFT_DIR)/passb4.f 

$(OBJ_DIR)/passb5.o :      $(FFT_DIR)/passb5.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passb5.o \
			   $(FFT_DIR)/passb5.f ${BOUNDS}
			    
$(OBJ_DIR)/passf.o :       $(FFT_DIR)/passf.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passf.o \
			   $(FFT_DIR)/passf.f ${BOUNDS}

$(OBJ_DIR)/passf2.o :      $(FFT_DIR)/passf2.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passf2.o \
			   $(FFT_DIR)/passf2.f 

$(OBJ_DIR)/passf3.o :      $(FFT_DIR)/passf3.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passf3.o \
			   $(FFT_DIR)/passf3.f ${BOUNDS}

$(OBJ_DIR)/passf4.o :      $(FFT_DIR)/passf4.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passf4.o \
			   $(FFT_DIR)/passf4.f 

$(OBJ_DIR)/passf5.o :      $(FFT_DIR)/passf5.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/passf5.o \
			   $(FFT_DIR)/passf5.f ${BOUNDS}

$(OBJ_DIR)/zfftb.o :       $(FFT_DIR)/zfftb.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/zfftb.o \
			   $(FFT_DIR)/zfftb.f 

$(OBJ_DIR)/zfftf.o :       $(FFT_DIR)/zfftf.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/zfftf.o \
			   $(FFT_DIR)/zfftf.f

$(OBJ_DIR)/zffti.o:        $(FFT_DIR)/zffti.f
			   $(COMPILE_FFT) -o $(OBJ_DIR)/zffti.o \
			   $(FFT_DIR)/zffti.f 
			   
				
### CQL3D_SETUP files:

$(OBJ_DIR)/basis_functions_m.o: $(CQL3D_SETUP_DIR)/basis_functions_m.f 
	                        $(COMPILE_CQL) -o $(OBJ_DIR)/basis_functions_m.o \
                                $(CQL3D_SETUP_DIR)/basis_functions_m.f ${BOUNDS}

$(OBJ_DIR)/f_expanded_m.o:    $(CQL3D_SETUP_DIR)/f_expanded_m.f 
	                      $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/f_expanded_m.o \
                              $(CQL3D_SETUP_DIR)/f_expanded_m.f ${BOUNDS}

$(OBJ_DIR)/global_data_m.o:   $(CQL3D_SETUP_DIR)/global_data_m.f 
	                      $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/global_data_m.o \
                              $(CQL3D_SETUP_DIR)/global_data_m.f ${BOUNDS}

$(OBJ_DIR)/CQL_kinds_m.o:     $(CQL3D_SETUP_DIR)/CQL_kinds_m.f 
	                      $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/CQL_kinds_m.o \
                              $(CQL3D_SETUP_DIR)/CQL_kinds_m.f ${BOUNDS}

$(OBJ_DIR)/vector_write_m.o:  $(CQL3D_SETUP_DIR)/vector_write_m.f 
	                      $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/vector_write_m.o \
                              $(CQL3D_SETUP_DIR)/vector_write_m.f ${BOUNDS}

$(OBJ_DIR)/read_cql3d.o:      $(CQL3D_SETUP_DIR)/read_cql3d.f 
	                      $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/read_cql3d.o \
                              $(CQL3D_SETUP_DIR)/read_cql3d.f ${BOUNDS}
			
#$(OBJ_DIR)/ceez.o:            $(CQL3D_SETUP_DIR)/ceez.f 
#	                      $(COMPILE_NOSAVE) -o $(OBJ_DIR)/ceez.o \
#                              $(CQL3D_SETUP_DIR)/ceez.f ${BOUNDS}
$(OBJ_DIR)/fitpack.o:            $(SRC_DIR)/fitpack.f 
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/fitpack.o \
                              $(SRC_DIR)/fitpack.f 


$(OBJ_DIR)/cubic_B_splines_v.o: $(CQL3D_SETUP_DIR)/cubic_B_splines_v.f 
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/cubic_B_splines_v.o \
                                $(CQL3D_SETUP_DIR)/cubic_B_splines_v.f ${BOUNDS}		
		
$(OBJ_DIR)/cql3d_setup.o:     $(CQL3D_SETUP_DIR)/cql3d_setup.f $(OBJ_DIR)/read_particle_f.o 
	                      $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/cql3d_setup.o \
                              $(CQL3D_SETUP_DIR)/cql3d_setup.f $(NETCDF) ${BOUNDS}
		
### DLG files:
		
$(OBJ_DIR)/read_particle_f.o:     $(DLG_DIR)/read_particle_f.f90 $(OBJ_DIR)/dlg.o $(OBJ_DIR)/eqdsk_dlg.o $(OBJ_DIR)/interp.o $(OBJ_DIR)/aorsa2din_mod.o
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/read_particle_f.o \
                              $(DLG_DIR)/read_particle_f.f90 $(NETCDF) ${BOUNDS}
ifeq ($(MACHINE),dlghp)
$(OBJ_DIR)/write_pf.o:     $(DLG_DIR)/write_pf.F ${OBJ_DIR}/write_pql.o
	                      $(COMPILE_DLG_R4) -DDLGHP -o $(OBJ_DIR)/write_pf.o \
                              $(DLG_DIR)/write_pf.F ${BOUNDS} ${PNETCDF} ${NETCDF}
else
$(OBJ_DIR)/write_pf.o:     $(DLG_DIR)/write_pf.F  ${OBJ_DIR}/write_pql.o
	                      $(COMPILE_DLG_R4) -DDLGHP -o $(OBJ_DIR)/write_pf.o \
                              $(DLG_DIR)/write_pf.F ${BOUNDS} ${PNETCDF} ${NETCDF}
endif
	
ifeq ($(MACHINE),dlghp)
$(OBJ_DIR)/write_pql.o:     $(DLG_DIR)/write_pql.F  
	                      $(COMPILE_DLG_R4) -DDLGHP -o $(OBJ_DIR)/write_pql.o \
                              $(DLG_DIR)/write_pql.F ${BOUNDS} ${PNETCDF}
else
$(OBJ_DIR)/write_pql.o:     $(DLG_DIR)/write_pql.F  
	                      $(COMPILE_DLG_R4) -DDLGHP -o $(OBJ_DIR)/write_pql.o \
                              $(DLG_DIR)/write_pql.F ${BOUNDS} ${PNETCDF}
endif
	
ifeq ($(MACHINE),dlghp)
$(OBJ_DIR)/gc_integrate.o:     $(DLG_DIR)/gc_integrate.F90 $(OBJ_DIR)/interp.o 
	                      $(COMPILE_DLG) -DUSE_DISLIN=1 -o $(OBJ_DIR)/gc_integrate.o \
                              $(DLG_DIR)/gc_integrate.F90 $(DISLIN) ${BOUNDS} 
else
$(OBJ_DIR)/gc_integrate.o:     $(DLG_DIR)/gc_integrate.F90 $(OBJ_DIR)/interp.o 
	                      $(COMPILE_DLG) -DUSE_DISLIN=0 -o $(OBJ_DIR)/gc_integrate.o \
                              $(DLG_DIR)/gc_integrate.F90 $(DISLIN) ${BOUNDS} 
endif
	
$(OBJ_DIR)/interp.o:     $(DLG_DIR)/interp.f90 $(OBJ_DIR)/gc_terms.o 
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/interp.o \
                              $(DLG_DIR)/interp.f90 ${BOUNDS}
	
$(OBJ_DIR)/gc_terms.o:     $(DLG_DIR)/gc_terms.f90 $(OBJ_DIR)/eqdsk_dlg.o $(OBJ_DIR)/constants.o 
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/gc_terms.o \
                              $(DLG_DIR)/gc_terms.f90 ${BOUNDS}
	
$(OBJ_DIR)/eqdsk_dlg.o:     $(DLG_DIR)/eqdsk_dlg.f90 $(OBJ_DIR)/dlg.o 
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/eqdsk_dlg.o \
                              $(DLG_DIR)/eqdsk_dlg.f90 ${BOUNDS}
	
$(OBJ_DIR)/dlg.o:     $(DLG_DIR)/dlg.f90 
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/dlg.o \
                              $(DLG_DIR)/dlg.f90 ${BOUNDS}
			
$(OBJ_DIR)/constants.o:     $(DLG_DIR)/constants.f90 
	                      $(COMPILE_DLG) -o $(OBJ_DIR)/constants.o \
                              $(DLG_DIR)/constants.f90 ${BOUNDS}

$(OBJ_DIR)/ranlib.o:        $(DLG_DIR)/ranlib.f 
	                     $(COMPILE) -o $(OBJ_DIR)/ranlib.o \
                             $(DLG_DIR)/ranlib.f ${BOUNDS}
	
make clean:
	rm $(EXEC) $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

