DEBUG ?= 1

ifeq ($(DEBUG), 1)
   CFLAGS=-Wall -g
   CxxFLAGS=-Wall -std=c++17 -g -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib
   CxxFLAGS_COMPLEX=-Wall -std=c++17 -g -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib
else
   CFLAGS=-Wall -O3 
   CxxFLAGS=-Wall -std=c++17 -O3 
endif


CCxx=mpicxx
CxxINCS=-I$(MFEM_DIR) -I$(MFEM_DIR)/linalg -I$(HYPRE_DIR)/src/hypre/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CxxLDIR=-L$(MFEM_DIR) -L$(HYPRE_DIR)/src/hypre/lib -L$(METIS_DIR) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/9
CxxLDIR_COMPLEX=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/9
CxxLIBS=-lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lptesmumps -lptscotchparmetisv3 -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lmfem -lHYPRE -lmetis -lrt -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd /usr/lib/x86_64-linux-gnu/liblapacke64.a /usr/lib/x86_64-linux-gnu/liblapack64.a -lpetscmat -lpetscsnes -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd -lgfortran -lc  -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp

CC=mpicc
CINCS=-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CLDIR=
CLIBS=-lm

fem.o: fem.cpp fem.hpp misc.hpp
	$(CCxx) $(CxxFLAGS) -c fem.cpp $(CxxINCS)

frequencyPlan.o: frequencyPlan.cpp frequencyPlan.hpp inputFrequency.h
	$(CCxx) $(CxxFLAGS) -c frequencyPlan.cpp $(CxxINCS)

jobrelated.o: jobrelated.cpp jobrelated.hpp 
	$(CCxx) $(CxxFLAGS) -c jobrelated.cpp $(CxxINCS)

keywordPair.o: keywordPair.cpp keywordPair.hpp misc.hpp
	$(CCxx) $(CxxFLAGS) -c keywordPair.cpp $(CxxINCS)

license.o: license.cpp license.hpp
	$(CCxx) $(CxxFLAGS) -c license.cpp $(CxxINCS)

mesh.o: mesh.cpp mesh.hpp misc.hpp
	$(CCxx) $(CxxFLAGS) -c mesh.cpp $(CxxINCS)

misc.o: misc.cpp misc.hpp
	$(CCxx) $(CxxFLAGS) -c misc.cpp $(CxxINCS)

OpenParEMmaterials.o: OpenParEMmaterials.cpp OpenParEMmaterials.hpp fem.hpp keywordPair.hpp misc.hpp
	$(CCxx) $(CxxFLAGS) -c OpenParEMmaterials.cpp $(CxxINCS)

path.o: path.cpp path.hpp keywordPair.hpp misc.hpp
	$(CCxx) $(CxxFLAGS) -c path.cpp $(CxxINCS)

petscErrorHandler.o: petscErrorHandler.cpp petscErrorHandler.hpp
	$(CCxx) $(CxxFLAGS) -c petscErrorHandler.cpp $(CxxINCS)

sourcefile.o: sourcefile.cpp sourcefile.hpp misc.hpp path.hpp
	$(CCxx) $(CxxFLAGS) -c sourcefile.cpp $(CxxINCS)

triplet.o: triplet.c triplet.h
	$(CC) $(CFLAGS) -c triplet.c $(CINCS)

libOpenParEMCommon.a: fem.o frequencyPlan.o jobrelated.o keywordPair.o license.o mesh.o misc.o OpenParEMmaterials.o path.o petscErrorHandler.o sourcefile.o triplet.o
	ar rcs libOpenParEMCommon.a fem.o frequencyPlan.o jobrelated.o keywordPair.o license.o mesh.o misc.o OpenParEMmaterials.o path.o petscErrorHandler.o sourcefile.o triplet.o

.PHONY: all clean install

all: libOpenParEMCommon.a

install:
	rm -f ../lib/libOpenParEMCommon.a
	cp libOpenParEMCommon.a ../lib/libOpenParEMCommon.a
	rm -f $(HOME)/lib/libOpenParEMCommon.a
	mv libOpenParEMCommon.a $(HOME)/lib/libOpenParEMCommon.a
	rm -f $(HOME)/bin/checkErrors.sh
	cp ../scripts/checkErrors.sh $(HOME)/bin/checkErrors.sh
	rm -f $(HOME)/listErrors
	cp ../scripts/listErrors $(HOME)/bin/listErrors
	rm -f $(HOME)/proj_search
	cp ../scripts/proj_search $(HOME)/proj_search

clean:
	rm -f libOpenParEMCommon.a
	rm -f fem.o
	rm -f frequencyPlan.o
	rm -f jobrelated.o
	rm -f keywordPair.o
	rm -f license.o
	rm -f mesh.o
	rm -f misc.o
	rm -f OpenParEMmaterials.o
	rm -f path.o
	rm -f petscErrorHandler.o
	rm -f sourcefile.o
	rm -f triplet.o
