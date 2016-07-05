# Makefile for C2Ray_3D.
#
# Author: Garrelt Mellema

# This Makefile can make different versions of C2Ray_3D.
# These versions differ in their parallelization and/or
# in their connection to specific N-body results.
#
# Note 1: Parallelization
# The parallelization intended is specified in the name
# of the executable: _omp means OpenMP (shared memory), 
# _mpi means MPI (distributed memory). Both can also be
# used at the same time (if your architecture supports
# it.
#
# Note 2: N-body module
# Different versions exist with different Nbody interfaces:
# pmfast - interface to older pmfast simulations
# cubep3m - interface to cubep3m simulations
# LG - interface to Local Group simulations (GADGET)
# Gadget - interface to LOFAR EoR GADGET simulations (not working)
#
# Note 3: Compiler & Flags
# The compiler is specified by the FC variable (MPIFC for the MPI
# compiler). We have only extensively used the Intel F90 compiler. 
# Support for other compilers will have to be added.
# Parts of the code need to know about the compiler, this is
# done through preprocessor statements. So when compiling with
# intel compiler, -DIFORT needs to be specified. Support for
# new compilers thus needs to be added in the code too.
#
# Note 4: Recompiling
# Some dependencies are through module parameters, and thus
# not recognized by make. Best practise is to run "make clean"
# before running "make".
#-------------------------------------------------------

# Compiler
FC = gfortran # GNU compiler
MPIFC = mpif90.openmpi # MPI compiler

# F90 options (gfortran)
GFORTFLAGS = -DGFORT -O3

# Processor dependent optimization
F90FLAGS1 = $(GFORTFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
OPENMP_FLAGS = -fopenmp -DMY_OPENMP # For Intel compiler

#-------------------------------------------------------
# Compiler
# Intel: best tested
#FC = ifort # Intel compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#IFORTFLAGS = -DIFORT -O0 -g -u -fpe0 -p
#IFORTFLAGS = -DIFORT -O3 -u -fpe0 -ipo -shared-intel #-check all -traceback
#IFORTFLAGS = -DIFORT -Ofast -fp-model fast=2 -u -fpe0 -ipo -shared-intel #-check all -traceback
#IFORTFLAGS = -DIFORT -Ofast -fp-model fast=2 -u -fpe0 #-check all -traceback
#IFORTFLAGS = -DIFORT -O3 -vec_report -u -fpe0 -ipo -mcmodel=medium -shared-intel #-check all -traceback

# Processor dependent optimization
#F90FLAGS1 = $(IFORTFLAGS) 
#F90FLAGS1 = -xCORE-AVX2 $(IFORTFLAGS) # Beskow
#F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xHost $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp -DMY_OPENMP # For Intel compiler
#NO_OPENMP_FLAGS = 

#-------------------------------------------------------
# Compiler
# Cray: 
#FC = ftn # Cray compiler
#MPIFC = ftn # Cray compiler

# F90 options (ifort)
#CRAYFLAGS = -O3 -hfp3 -rm -DCRAY #-check all -traceback
# Processor dependent optimization
#F90FLAGS1 = $(CRAYFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI # For the code
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -DMY_OPENMP # For the code
#NO_OPENMP_FLAGS = -Othread0
#-------------------------------------------------------

# Compiler
# Sun: problems with constant definition. Cannot have sqrt in constant
# definition.
#FC = f95 # Sun compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#SUNFLAGS = -O3 -DSUN
# Processor dependent optimization
#F90FLAGS1 = $(SUNFLAGS) 
#F90FLAGS1 = -xW $(SUNFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp -DMY_OPENMP # For Sun compiler

#-------------------------------------------------------

# PGI compiler
#FC = pf90
#MPIFC = mpif77
#MPIFC = mpif90

# F90 options (pgi)
#PGIFLAGS = -O3 -fast -DPGI
#F90FLAGS1 = -tp barcelona-64  $(PGIFLAGS) # ranger processors

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -mp -DMY_OPENMP

#-------------------------------------------------------

# Other F90 compilers
#FC = f90

#-------------------------------------------------------

#LDR     = $(F90)

LDFLAGS = $(F90FLAGS) #-L/afs/astro.su.se/pkg/intel/Compiler/11.1/056/lib/intel64/
LIBS = #-lirc -limf

#-------------------------------------------------------

UTILS=mrgrnk.o ctrper.o romberg.o report_memory.o

CONSTANTS = mathconstants.o cgsconstants.o  cgsphotoconstants.o  cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o

#--------TEST----------------------------------------------------------------

C2Ray_3D_test_periodic: F90=$(FC)
C2Ray_3D_test_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_test_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o mpi.o clocks.o file_admin.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_test_periodic_omp: F90=$(FC)
C2Ray_3D_test_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_test_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

#--------PMFAST---------------------------------------------------------------

C2Ray_3D_pmfast_periodic: F90=$(FC)
C2Ray_3D_pmfast_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_pmfast_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic_compr: F90=$(FC)
C2Ray_3D_pmfast_periodic_compr: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_pmfast_periodic_compr: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o pmfast.o grid.o tped.o mat_ini_pmfast_compr.o sourceprops_pmfast_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast_compr.o sourceprops_pmfast_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o

C2Ray_3D_pmfast_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_pmfast_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_pmfast_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o mpi.o clocks.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_pmfast_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_pmfast_periodic_omp_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o mpi.o clocks.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic_omp: F90=$(FC)
C2Ray_3D_pmfast_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_pmfast_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

#--------CUBEP3M-------------------------------------------------------------

C2Ray_3D_cubep3m_periodic: F90=$(FC)
C2Ray_3D_cubep3m_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_cubep3m_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o src_subgrid.o mergesrc.o cooling.o radiation.o jLWgreen.o cosmology.o time_ini.o getjLW.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o src_subgrid.o mergesrc.o cooling.o radiation.o jLWgreen.o cosmology.o getjLW.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS) $(NO_OPENMP_FLAGS)
C2Ray_3D_cubep3m_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(CONSTANTS) $(UTILS) sizes.o mpi.o clocks.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_omp: F90=$(FC)
C2Ray_3D_cubep3m_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) 
C2Ray_3D_cubep3m_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(CONSTANTS) $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_omp_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(CONSTANTS) $(UTILS) sizes.o mpi.o clocks.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_compr_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_compr_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_compr_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m_compr.o sourceprops_cubep3m_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o mpi.o clocks.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m_compr.o sourceprops_cubep3m_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o

#--------LG------------------------------------------------------------------
 
C2Ray_3D_LG_periodic: F90=$(FC)
C2Ray_3D_LG_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_LG_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_LG_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_LG_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_LG_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o mpi.o clocks.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_LG_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_LG_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_LG_periodic_omp_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o clocks.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o mpi.o clocks.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

C2Ray_3D_LG_periodic_omp: F90=$(FC)
C2Ray_3D_LG_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) 
C2Ray_3D_LG_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray.o

#--------GADGET-------------------------------------------------------------------------------

C2Ray_3D_Gadget_periodic: F90=$(FC)
C2Ray_3D_Gadget_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_Gadget_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray_GadgetTest.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray_GadgetTest.o

C2Ray_3D_Gadget_periodic_omp: F90=$(FC)
C2Ray_3D_Gadget_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_Gadget_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o clocks.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray_GadgetTest.o
	$(F90) $(F90FLAGS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o clocks.o file_admin.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve8.o output.o C2Ray_GadgetTest.o

clean : 
	rm -f *.o *.mod *.l *.il

.f.o:
	$(F90) -c $(F90FLAGS) $<

.f90.o:
	$(F90) -c $(F90FLAGS) $<

.F90.o:
	$(F90) -c $(F90FLAGS) $<

f.mod:
	$(F90) -c $(F90FLAGS) $<

.SUFFIXES: .f90 .F90 .mod .o


