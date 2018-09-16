#### PRESUMED STRUCTURE
#		>>build        !(You are here to call the "make" command)
#		>>randomField  !(This is the git repository)
#			-BLAS
#			-SRC_EXE
#			-SRC_LIB
#

#### PATH TO HDF5 AND MPI LIBRARY AND INCLUDES - to be modified by user (to be completed)

#LIBHDF5     = -L/gpfs/opt/libraries/phdf5/1.8.17-IntelMPI-bis/lib/ -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran
#INCLUDEHDF5 = -I/gpfs/opt/libraries/phdf5/1.8.17-IntelMPI-bis/include
#
#
#LIBMPI      = -L/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64
#INCLUDEMPI  = -I/gpfs/opt/compilers/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include
#
#LIBFFTW     = -L/gpfs/opt/libraries/fftw/3.3.4/lib -lfftw3_mpi -lfftw3
#INCLUDEFFTW = -I/gpfs/opt/libraries/fftw/3.3.4/include

LIBHDF5     = -L/home/carvalho/LOCAL_LIB/lib/ -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran
INCLUDEHDF5 = -I/home/carvalho/LOCAL_LIB/include
LIBFFTW     = -L/home/carvalho/LOCAL_LIB/lib -lfftw3_mpi -lfftw3
INCLUDEFFTW = -I/home/carvalho/LOCAL_LIB/include

LIBMPI      = -L/gpfs/opt/compilers/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/lib -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64
INCLUDEMPI  = -I/gpfs/opt/compilers/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/include

EXEC2 = statistics.exe
FC   = mpiifort
#FFLAGS = -g -cpp
FFLAGS = -xCORE-AVX2 -O3 -no-prec-div

#### FROM THIS POINT ON THE USER SHOULD NOT MODIFY
RM := rm -rf

#Specifying the sources "f90" and "f"
SRCS = $(wildcard *.f90 *.f) 

#Create a ".f90" for each source
#OBJS = $(SRCS:.f90=.o) #SYNTAX NOT WORKING

STAT_OBJ += ./main_Stat.o \
./statistics_RF.o \
./type_STAT.o \
./charFunctions.o \
./constants_RF.o \
./common_variables_RF.o \
./displayCarvalhol.o \
./hdf5_RF.o \
./write_Log_File.o \
./readFile_RF.o \
./math_RF.o \
./dgemm.o \
./lsame.o \
./xerbla.o \
./dger.o \
./systemUt_RF.o \
./fftw3.o

#./type_TEST.o \
#./test_func_RF.o \

LIBS = $(LIBHDF5) $(LIBMPI) $(LIBFFTW)
INCLUDE = $(INCLUDEHDF5) $(INCLUDEMPI) $(INCLUDEFFTW)

#Dependencies
main_Stat.o          : charFunctions.o \
                       common_variables_RF.o \
		       constants_RF.o \
		       displayCarvalhol.o \
	               hdf5_RF.o \
		       math_RF.o \
                       write_Log_File.o \
                       readFile_RF.o \
                       systemUt_RF.o \
		       statistics_RF.o \
		       type_STAT.o	           				  
hdf5_RF.o            : displayCarvalhol.o \
                       math_RF.o \
                       write_Log_File.o
statistics_RF.o      : displayCarvalhol.o \
		       math_RF.o \
		       write_Log_File.o \
		       type_STAT.o
math_RF.o            : displayCarvalhol.o \
		       write_Log_File.o \
		       constants_RF.o
systemUt_RF.o        : write_Log_File.o
write_Log_File.o     : charFunctions.o
charFunctions.o      : constants_RF.o
common_variables_RF.o      : constants_RF.o
displayCarvalhol.o   : constants_RF.o \
                       write_Log_File.o
dgemm.o              : lsame.o xerbla.o
dger.o               : xerbla.o
# Making all the ".o" from the ".f90"
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../BLAS/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '
	
# All Target
all: statistics
#all: randomField

# Tool invocations
	
statistics: $(STAT_OBJ)
		@echo 'Building target: $@'
		@echo 'Invoking: Fortran Linker'
		$(FC) -o $(EXEC2) $(FFLAGS) $(STAT_OBJ) $(INCLUDE) $(LIBS)
		@echo 'Finished building target: $@'
		@echo ' '

# Other Targets
clean:
		-$(RM) $(EXEC2) $(STAT_OBJ)  *.mod
		-@echo 'Finished excluding'

.PHONY: all clean dependents
.SECONDARY:
