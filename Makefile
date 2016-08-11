TARGET=Cc4s
INSTALL=~/bin/cc4s
include config.mk
VERSION:=$(shell git describe --all --dirty --long)
DATE:=$(shell git log -1 --format="%cd")
COMPILER_VERSION:=$(shell ${CXX} --version | head -n 1)
# location of the Cyclops Tensor Framework library
CTF=../ctf
OPTIMIZE=${CXXOPTIMIZE} -O3
COPTIONS=${CXXOPTIONS} -std=c++0x -Wall -fmax-errors=3 -g \
-D_POSIX_C_SOURCE=200112L \
-D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1 -DCC4S_VERSION=\"${VERSION}\" \
"-DCC4S_DATE=\"${DATE}\"" \
"-DCOMPILER_VERSION=\"${COMPILER_VERSION}\""

# primary target
cc4s: bin/${TARGET}

doc:
	doxygen

install: bin/${TARGET}
	mkdir -p ${INSTALL}
	cp bin/${TARGET} ${INSTALL}

OBJECTS= \
obj/Options.o \
obj/util/Log.o \
obj/util/Timer.o \
obj/util/FlopsCounter.o \
obj/util/DryTensor.o \
obj/math/MathFunctions.o \
obj/math/ComplexTensor.o \
obj/math/RandomTensor.o \
obj/math/IterativePseudoInverse.o \
obj/math/CanonicalPolyadicDecomposition.o \
obj/math/RegularizedAlternatingLeastSquares.o \
obj/PerturbationTensor.o \
obj/Chi.o obj/CoulombIntegrals.o obj/Amplitudes.o \
obj/BinaryFtodReader.o obj/TextFtodReader.o \
obj/Data.o \
obj/Parser.o \
obj/mixers/Mixer.o \
obj/mixers/LinearMixer.o \
obj/algorithms/Algorithm.o \
obj/algorithms/ClusterDoublesAlgorithm.o \
obj/algorithms/ClusterSinglesDoublesAlgorithm.o \
obj/algorithms/ParticleHoleCoulombIntegrals.o \
obj/algorithms/ParticleHoleCoulombVertexReader.o \
obj/algorithms/DrccdEnergyFromCoulombIntegrals.o \
obj/algorithms/DrccdEnergyFromCoulombVertex.o \
obj/algorithms/SliceCoulombVertex.o \
obj/algorithms/CoulombVertexDecomposition.o \
obj/algorithms/ParticleHoleCoulombVertexDecomposition.o \
obj/algorithms/GenerateRandomMatrix.o \
obj/algorithms/GenerateRandomComplexMatrix.o \
obj/algorithms/GenerateRandomTensor.o \
obj/algorithms/TensorReader.o \
obj/algorithms/TensorWriter.o \
obj/algorithms/TensorContraction.o \
obj/algorithms/TensorSum.o \
obj/algorithms/TensorNorm.o \
obj/algorithms/ComplexTensorWriter.o \
obj/algorithms/ComplexTensorContraction.o \
obj/algorithms/ComplexTensorSum.o \
obj/algorithms/ComplexTensorNorm.o \
obj/algorithms/FromComplexTensor.o \
obj/algorithms/PartitionTensor.o \
obj/algorithms/CoulombVertexReader.o \
obj/algorithms/CoulombIntegralsFromVertex.o \
obj/algorithms/Mp2EnergyFromCoulombIntegrals.o \
obj/algorithms/CcdEnergyFromCoulombIntegrals.o \
obj/algorithms/DcdEnergyFromCoulombIntegrals.o \
obj/algorithms/CcdEnergyFromCoulombFactors.o \
obj/algorithms/DcdEnergyFromCoulombFactors.o \
obj/algorithms/CcsdEnergyFromCoulombIntegrals.o \
obj/algorithms/DcsdEnergyFromCoulombIntegrals.o \
obj/algorithms/Mp2EnergyMatrixFromCoulombIntegrals.o \
obj/algorithms/EnergyMatrixFromDoublesAmplitudes.o \
obj/algorithms/ReduceEnergyMatrix.o \
obj/algorithms/ReduceCoulombVertex.o \
obj/algorithms/ReduceParticleHoleCoulombVertex.o \
obj/algorithms/CcsdEnergyFromCoulombFactors.o \
obj/algorithms/DcsdEnergyFromCoulombFactors.o \
obj/algorithms/CcdEnergyFromSlicedIntegrals.o \
obj/algorithms/DcdEnergyFromSlicedIntegrals.o \
obj/algorithms/CcsdEnergyFromSlicedIntegrals.o \
obj/algorithms/DcsdEnergyFromSlicedIntegrals.o \
obj/algorithms/SingleParticleOccupancies.o \
obj/algorithms/PolarizabilityFromCoulombVertex.o \
# dependencies
obj/${TARGET}.o: ${OBJECTS}

clean:
	rm -rf bin/*
	rm -rf obj/*

# compile object files
obj/%.o: src/%.cxx src/%.hpp
	mkdir -p $(dir $@)
	${CXX} ${COPTIONS} ${OPTIMIZE} -c src/$*.cxx -o $@ -I${CTF}/include -Isrc

# compile and link executable
bin/%: obj/%.o
	${CXX} ${COPTIONS} ${OPTIMIZE} ${OBJECTS} obj/${TARGET}.o -o $@ -I${CTF}/include/ -L${CTF}/lib -lctf ${LIBS}

