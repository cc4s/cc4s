TARGET=Cc4s
INSTALL=~/bin/cc4s
VERSION:=$(shell git describe --all --dirty --long)
DATE:=$(shell git log -1 --format="%cd")
# location of the Cyclops Tensor Framework library
CTF=../ctf
include config.mk
OPTIMIZE=${CXXOPTIMIZE} -O3
COPTIONS=${CXXOPTIONS} -std=c++0x -Wall -fmax-errors=3 \
-D_POSIX_C_SOURCE=200112L \
-D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1 -DCC4S_VERSION=\"${VERSION}\" \
"-DCC4S_DATE=\"${DATE}\""
#LIBS=-lblas -lgfortran

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
obj/util/MathFunctions.o \
obj/util/ComplexTensor.o \
obj/util/RandomTensor.o \
obj/util/IterativePseudoInverse.o \
obj/PerturbationTensor.o \
obj/Data.o \
obj/Algorithm.o \
obj/Parser.o \
obj/RalsFtodRankDecomposition.o \
obj/Chi.o obj/CoulombIntegrals.o obj/Amplitudes.o \
obj/BinaryFtodReader.o obj/TextFtodReader.o \
obj/ParticleHoleCoulombIntegrals.o \
obj/ParticleHoleCoulombVertexReader.o \
obj/Mp2EnergyFromCoulombIntegrals.o \
obj/DrccdEnergyFromCoulombIntegrals.o \
obj/DrccdEnergyFromCoulombVertex.o \
obj/RalsParticleHoleCoulombVertexDecomposition.o \
obj/CtfTest.o \
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

