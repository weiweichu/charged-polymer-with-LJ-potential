optstd=c++98
optstd=c++0x# for "long long" type support, among others
src_dir=.
cmp_def=-DGREEN_SYMMETRY
cmp_def=-DGREEN_CONSISTENCE
cmp_def=-DGREEN_NEW
cxxflags= -O3 -ffast-math -Wall -Winline -std=$(optstd) -pedantic -DHAVE_INLINE $(cmp_def) -lgsl -lgslcblas -lfftw3 

cxx=g++ $(cxxflags) -I$(src_dir)
exe=gmc

#cxx=mpicxx $(cxxflags) -I$(src_dir) -DUTIL_MPI
#exe=gmc_p

constants=util/Constants
log=util/Log
exception=util/Exception
vector=util/Vector
intvector=util/IntVector
mersenne=util/Mersenne
filemaster=util/FileMaster
mpitraits=util/mpi/MpiTraits
mpistructbuilder=util/mpi/MpiStructBuilder
util_obj=$(constants).o $(log).o $(exception).o $(vector).o $(intvector).o $(mersenne).o $(filemaster).o $(mpitraits).o $(mpistructbuilder).o

grid=grid/Grid
gridmass=grid/GridMass
gridmasscharge=grid/GridMassCharge
poissongreenfunc=grid/PoissonGreenFunc
grid_obj=$(grid).o $(gridmass).o $(gridmasscharge).o $(poissongreenfunc).o

ewald=charge/Ewald
charge_obj=$(ewald).o

move=move/Move
beadmove=move/BeadMove
volumemove=move/VolumeMove
chainflipmove=move/ChainFlipMove
reptationmove=move/ReptationMove
semigrandsaltmove=move/SemigrandSaltMove
move_obj=$(move).o $(beadmove).o $(volumemove).o $(chainflipmove).o $(reptationmove).o $(semigrandsaltmove).o

oprm=oprm/OrderParameter
lamoprm=oprm/LamOrderParameter
hexoprm=oprm/HexOrderParameter
oprm_obj=$(oprm).o $(lamoprm).o $(hexoprm).o

diag=diag/Diagnosis
structurefactor=diag/StructureFactor
solventchemicalpotential=diag/SolventChemicalPotential
saltchemicalpotential=diag/SaltChemicalPotential
chainchemicalpotential=diag/ChainChemicalPotential
chainchemicalpotentialinsertion=diag/ChainChemicalPotentialInsertion
ionnumberfraction=diag/IonNumberFraction
diag_obj=$(diag).o $(structurefactor).o $(solventchemicalpotential).o $(saltchemicalpotential).o $(ionnumberfraction).o $(chainchemicalpotential).o $(chainchemicalpotentialinsertion).o

particle=simulation/Particle
histogram=simulation/Histogram
histogramweight=simulation/HistogramWeight
system=simulation/System
mcsystem=simulation/McSystem
simulation=simulation/Simulation
simulation_obj=$(particle).o $(histogram).o $(histogramweight).o $(system).o $(mcsystem).o $(simulation).o

obj=$(util_obj) $(grid_obj) $(charge_obj) $(move_obj) $(oprm_obj) $(diag_obj) $(simulation_obj)

.PHONY: clean cleanall

clean:
	-rm -f $(obj)

cleanall: clean
	-rm -f $(exe)

all: $(obj) $(exe)

$(exe): gmc.cpp $(obj)
	$(cxx) -o $(exe) gmc.cpp $(obj)

%.o: %.cpp %.h
	$(cxx) -o $@ -c $<

