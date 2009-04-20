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
# The Makefile does NOT automatically set the right flags
# for the executable you want to make. 
# If you want to compile for OpenMP, you have to add the
# the OpenMP flags (OPENMP_FLAGS) to F90FLAGS manually 
# before running make.
# If you want to compile for MPI you have to add the MPI
# flags (MPI_FLAGS) to F90FLAGS, and change the compiler (F90) 
# manually before running make.
#
# Note 2: N-body module
# Different versions exist with different Nbody interfaces:
# pmfast - interface to older pmfast simulations
# cubep3m - interface to cubep3m simulations
# LG - interface to Local Group simulations (GADGET)
# Gadget - interface to LOFAR EoR GADGET simulations (not working)
#
# Note 3: Compiler & Flags
# The compiler is specified by the F90 variable. We have only
# extensively used the Intel F90 compiler. Support for other
# compilers will have to added.
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
FC = ifort # Intel compiler
MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
F90FLAGS1 = $(IFORTFLAGS) 
#F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xO $(IFORTFLAGS) 
F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
MPI_FLAGS = -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
OPENMP_FLAGS = -openmp # For Intel compiler

#-------------------------------------------------------

# PGI compiler
#F90 = pf90
#F90 = mpif77
#F90 = mpif90

# F90 options (pgi)
#PGIFLAGS = -O3 -fast -DPGI
#F90FLAGS1 = -tp barcelona-64  $(PGIFLAGS) # ranger processors

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -mp 

#-------------------------------------------------------

# Other F90 compilers
#F90 = f90

#F90FLAGS = -O3 -u -fpe0 -vec_report -ipo #Lobster 6
#F90FLAGS = -O4 -omp -fpe1# -fpe4 Compaq
#F90FLAGS  = -O3 -openmp

#F90FLAGS = -g -O0 -openmp #seg fault
#F90FLAGS = -g -O0 
#F90FLAGS = -u -fpe0
#F90FLAGS = -O3 -u -ipo -fpe0
#F90FLAGS = -O3 -ipo
#F90FLAGS = -fpe0 -g -u 
#F90FLAGS = -O3 -vec_report -u -ipo
#F90FLAGS = -xW -O3 -openmp -vec_report -u -ipo -fpe0
#F90FLAGS =  -Wv,-Pprocs,4 -O3 -cpu:opteron

#-------------------------------------------------------

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS)
LIBS = 

#-------------------------------------------------------

UTILS=mrgrnk.o ctrper.o romberg.o

CONSTANTS = mathconstants.o cgsconstants.o  cgsphotoconstants.o  cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o

C2Ray_3D_test_periodic: F90=$(FC)
C2Ray_3D_test_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_test_periodic_box: F90=$(FC)
C2Ray_3D_test_periodic_box: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test_periodic_box: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve5.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o test.o grid.o tped.o mat_ini_test.o sourceprops_test.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve5.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic: F90=$(FC)
C2Ray_3D_pmfast_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_pmfast_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic_compr: F90=$(FC)
C2Ray_3D_pmfast_periodic_compr: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_pmfast_periodic_compr: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o pmfast.o grid.o tped.o mat_ini_pmfast_compr.o sourceprops_pmfast_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast_compr.o sourceprops_pmfast_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o

C2Ray_3D_pmfast_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_pmfast_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_pmfast_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_pmfast_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_pmfast_periodic_omp_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o

C2Ray_3D_pmfast_periodic_omp: F90=$(FC)
C2Ray_3D_pmfast_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_pmfast_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic: F90=$(FC)
C2Ray_3D_cubep3m_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_cubep3m_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_box_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_box_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_box_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve5.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve5.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_nbox_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_nbox_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_nbox_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve6.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve6.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_compr_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_compr_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_compr_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m_compr.o sourceprops_cubep3m_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m_compr.o sourceprops_cubep3m_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o

C2Ray_3D_cubep3m_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_omp_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o

C2Ray_3D_cubep3m_periodic_omp: F90=$(FC)
C2Ray_3D_cubep3m_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) 
C2Ray_3D_cubep3m_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o

C2Ray_3D_LG_periodic: F90=$(FC)
C2Ray_3D_LG_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_LG_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_LG_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_LG_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_LG_periodic_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray.o

C2Ray_3D_LG_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_LG_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_LG_periodic_omp_mpi: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o mpi.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o

C2Ray_3D_LG_periodic_omp: F90=$(FC)
C2Ray_3D_LG_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) 
C2Ray_3D_LG_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray.o

C2Ray_3D_Gadget_periodic: F90=$(FC)
C2Ray_3D_Gadget_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_Gadget_periodic: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray_GadgetTest.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_periodic.o output.o C2Ray_GadgetTest.o

C2Ray_3D_Gadget_periodic_omp: F90=$(FC)
C2Ray_3D_Gadget_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_Gadget_periodic_omp: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o no_mpi.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray_GadgetTest.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) sizes.o no_mpi.o file_admin.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics.o evolve4_omp_periodic.o output.o C2Ray_GadgetTest.o

clean : 
	rm -f *.o *.mod *.l *.il

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o


