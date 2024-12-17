
SHELL=/bin/sh
.SUFFIXES:	
.SUFFIXES:	.F .o
.SUFFIXES:  .F90 .o
OBJS=	extract_apar_wrf.o \
	module_cfradial_output.o \
	kwm_date_utilities.o \
	module_llxy.o \
	module_external_attitude.o \
	module_access_wrf.o \
	module_crsim_wrapper.o \
	module_quicksort.o \
	crsim_mod.o \
	phys_param_mod.o \
	module_scanning.o \
	module_aircraft.o \
	module_geometry.o \
	wrf_var_mod.o \
	crsim_luts_mod.o \
	ReadConfParameters.o  \
	module_configuration.o \
	module_timing_utilities.o

#####  GNU  ###################################################################
LD = gfortran
FC = gfortran
FORTRAN_FREEFORM = -ffree-form -ffree-line-length-none
NETCDF_INC = $(NETCDF)/include
NETCDF_LIB = $(NETCDF)/lib
FFLAGS =  -I$(NETCDF_INC) $(FORTRAN_FREEFORM) -O3
LDLIBS =	-L$(NETCDF_LIB) -lnetcdff
DEBUG  = 	# -g -C -fbounds-check -O0

#####  GNU + MPI  #############################################################
# LD = mpif90 
# FC = mpif90 -D_PARALLEL_
# FORTRAN_FREEFORM = -ffree-form -ffree-line-length-none
# NETCDF_INC = $(NETCDF)/include
# NETCDF_LIB = $(NETCDF)/lib
# FFLAGS = -I$(NETCDF_INC) $(FORTRAN_FREEFORM) -O3
# LDLIBS = -L$(NETCDF_LIB) -lnetcdff
# DEBUG  = # -g -C -fbounds-check -O0

#####  Intel  #################################################################
# LD     = 	ifort
# FC     = 	ifort
# FORTRAN_FREEFORM = -free
# NETCDF_INC = $(NETCDF)/include
# NETCDF_LIB = $(NETCDF)/lib
# FFLAGS =	-I$(NETCDF_INC) $(FORTRAN_FREEFORM) -O3 # -r8
# LDLIBS =	-L$(NETCDF_LIB) -lnetcdff
# DEBUG  = 	# -g -C -O0

#####  Intel + MPI  ###########################################################
# LD = mpif90
# FC = mpif90 -D_PARALLEL_
# FORTRAN_FREEFORM = -free
# NETCDF_INC = $(NETCDF)/include
# NETCDF_LIB = $(NETCDF)/lib
# FFLAGS =	-I$(NETCDF_INC) $(FORTRAN_FREEFORM) -O3 # -r8
# LDLIBS =	-L$(NETCDF_LIB) -lnetcdff
# DEBUG  = 	# -g -C -O0

#####  nvfortran  #############################################################
# LD = nvfortran
# FC = nvfortran
# FORTRAN_FREEFORM = -Mfree
# NETCDF_INC = $(NETCDF)/include
# NETCDF_LIB = $(NETCDF)/lib
# FFLAGS = $(FORTRAN_FREEFORM) -O3
# LDLIBS =	-L$(NETCDF_LIB) -lnetcdff
# DEBUG = # -g -C -O0

#####  nvfortran + MPI  #######################################################
# LD = mpif90 
# FC = mpif90 -D_PARALLEL_
# FORTRAN_FREEFORM = -Mfree
# NETCDF_INC = $(NETCDF)/include
# NETCDF_LIB = $(NETCDF)/lib
# FFLAGS = $(FORTRAN_FREEFORM) -O3
# LDLIBS =	-L$(NETCDF_LIB) -lnetcdff
# DEBUG = # -g -C -O0

#####  Using nf-config from installation of NetCDF  ###########################
# LD = $(shell nf-config --fc)
# FC = $(shell nf-config --fc)
# ifeq ($(FC),ifort)
# 	FFLAGS = -free -O3 $(shell nf-config --fflags)
# endif
# ifeq ($(FC),gfortran)
# 	FFLAGS = -ffree-form -ffree-line-length-none -O3 $(shell nf-config --fflags)
# endif
# LDLIBS = $(shell nf-config --flibs)

###############################################################################

RM = 	rm -f
CMD=	aospre

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build $(CMD).
#
all:	$(CMD)

aospre:	$(OBJS)
	$(LD) -o $(@) $(OBJS) $(LDLIBS)

.f.o:
	$(FC) -c $(FFLAGS) $(DEBUG) $*.f

.c.o:
	$(CC) -c $(CCFLAGS) $*.c

.F.o:
	$(FC) -c $(FFLAGS) $(DEBUG) $(FCPPFLAGS) $*.F

.F90.o:
	$(FC) -c $(FFLAGS) $(DEBUG) $(FCPPFLAGS) $*.F90

clean:
	$(RM) *.o *.mod *~ $(CMD)
#

extract_apar_wrf.o : module_cfradial_output.o
extract_apar_wrf.o : module_access_wrf.o
extract_apar_wrf.o : module_scanning.o
extract_apar_wrf.o : module_aircraft.o
extract_apar_wrf.o : module_llxy.o
extract_apar_wrf.o : module_geometry.o
extract_apar_wrf.o : module_external_attitude.o
extract_apar_wrf.o : module_crsim_wrapper.o
extract_apar_wrf.o : crsim_mod.o
extract_apar_wrf.o : module_configuration.o
extract_apar_wrf.o : module_timing_utilities.o
module_access_wrf.o: kwm_date_utilities.o
module_access_wrf.o: module_configuration.o
module_access_wrf.o: module_llxy.o
module_access_wrf.o: module_scanning.o
module_cfradial_output.o : kwm_date_utilities.o
module_cfradial_output.o : module_llxy.o
module_cfradial_output.o : module_configuration.o
module_aircraft.o: module_llxy.o
module_aircraft.o: module_timing_utilities.o
module_aircraft.o: module_external_attitude.o
module_aircraft.o: module_configuration.o
module_aircraft.o: module_access_wrf.o
module_geometry.o: module_llxy.o
module_geometry.o: module_configuration.o
module_scanning.o: module_configuration.o
module_llxy.o: module_configuration.o
module_external_attitude.o: module_configuration.o
module_crsim_wrapper.o: crsim_mod.o
module_crsim_wrapper.o: phys_param_mod.o
module_crsim_wrapper.o: wrf_var_mod.o
module_crsim_wrapper.o: crsim_luts_mod.o
module_crsim_wrapper.o: module_cfradial_output.o
module_crsim_wrapper.o: module_configuration.o
module_crsim_wrapper.o: module_quicksort.o
module_quicksort.o: module_configuration.o
module_quicksort.o: quicksort_inline.inc
module_quicksort.o: quicksort_inline_index.inc
crsim_mod.o: phys_param_mod.o
ReadConfParameters.o: crsim_mod.o
