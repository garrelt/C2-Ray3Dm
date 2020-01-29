
# Compiler
FC = ifort # Intel compiler
MPIFC = mpif90 # MPI compiler

# FFTW Chosun Teragon
FFTW3 =/usr/local/fftw-3.2.2

# FFTW Ranger local
#FFTW3 =/share/home/01113/tg803913/fftw_3.2.2

FFTW3_INC =$(FFTW3)/include
FFTW3_LIB =$(FFTW3)/lib
# place FFTW3_FLAG after .o files when linking fails.
FFTW3_FLAG = -I$(FFTW3_INC) -L$(FFTW3_LIB) -lfftw3

# F90 options (ifort)
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT -shared-intel #-check all -traceback
IFORTFLAGS = -O3 -vec_report -u -fpp -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -DIFORT -check all -traceback
#IFORTFLAGS = -O2 -fpp  # Chosun Teragon
# Processor dependent optimization
F90FLAGS1 = $(IFORTFLAGS) 
F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xO $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
MPI_FLAGS = -DMPI # 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#MPI_FLAGS = -DMPI -I/engrid/enhpc/mpich/intel/include
OPENMP_FLAGS = -openmp -D_OPENMP# For Intel compiler

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
#OPENMP_FLAGS = -mp 

#-------------------------------------------------------

# Other F90 compilers
#FC = f90
F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)                 # MPI only

#-------------------------------------------------------

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS)
LIBS = 

#-------------------------------------------------------

UTILS= spline_simple_module.o

CONSTANTS = mathconstants.o cgsconstants.o cgsastroconstants.o cosmoparms.o


morph_dens: F90=$(FC)
morph_dens: F90FLAGS = $(F90FLAGS1)
morph_dens: precision.o $(CONSTANTS) $(UTILS) cosmology.o vtk_2d.o morph_dens.o
	$(F90) precision.o $(UTILS) cosmology.o vtk_2d.o morph_dens.o $(OPTIONS) -o $@ 

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


