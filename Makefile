TARGET=Cc4s
VERSION:=$(shell git describe --all --dirty --long)
DATE:=$(shell git log -1 --format="%cd")
# location of the Cyclops Tensor Framework library
CTF=../ctf
#TODO: use configuration files
CXX=mpicxx
OPTIMIZE=-O3
COPTIONS=-std=c++0x -fopenmp -Wall -g -fmax-errors=3 -D_POSIX_C_SOURCE=200112L \
-D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1 -DCC4S_VERSION=\"${VERSION}\" \
"-DCC4S_DATE=\"${DATE}\""
LIBS=-lblas -lgfortran

# primary target
cc4s: bin/${TARGET}

OBJECTS= \
obj/Options.o \
obj/util/Log.o \
obj/util/MathFunctions.o \
obj/util/ComplexTensor.o \
obj/util/RandomTensor.o \
obj/util/IterativePseudoInverse.o \
obj/PerturbationTensor.o \
obj/Algorithm.o \
obj/FtodRankDecomposition.o \
obj/RalsFtodRankDecomposition.o \
obj/CrossEntropyFtodRankDecomposition.o \
obj/Chi.o obj/CoulombIntegrals.o obj/Amplitudes.o \
obj/BinaryFtodReader.o obj/TextFtodReader.o obj/${TARGET}.o \
obj/ParticleHoleCoulomb.o \
obj/ParticleHoleCoulombVertexReader.o
# dependencies
obj/${TARGET}.o: \
obj/Options.o \
obj/util/Log.o \
obj/util/MathFunctions.o \
obj/util/ComplexTensor.o \
obj/util/RandomTensor.o \
obj/util/IterativePseudoInverse.o \
obj/PerturbationTensor.o \
obj/Algorithm.o \
obj/FtodRankDecomposition.o \
obj/RalsFtodRankDecomposition.o \
obj/CrossEntropyFtodRankDecomposition.o \
obj/Chi.o obj/CoulombIntegrals.o obj/Amplitudes.o \
obj/BinaryFtodReader.o obj/TextFtodReader.o \
obj/ParticleHoleCoulomb.o \
obj/ParticleHoleCoulombVertexReader.o

clean:
	rm -rf bin/*
	rm -rf obj/*

# compile object files
obj/%.o: src/%.cxx src/%.hpp
	mkdir -p $(dir $@)
	${CXX} ${COPTIONS} ${OPTIMIZE} -c src/$*.cxx -o $@ -I${CTF}/include -Isrc

# compile and link executable
bin/%: obj/%.o
	${CXX} ${COPTIONS} ${OPTIMIZE} ${OBJECTS} -o $@ -I${CTF}/include/ -L${CTF}/lib -lctf ${LIBS}

