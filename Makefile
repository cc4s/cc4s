# location of the Cyclops Tensor Framework library
CTF=../ctf
#TODO: use configuration files
CXX=mpicxx
OPTIMIZE=-O3
COPTIONS=-std=c++0x -fopenmp -Wall -g -fmax-errors=3 -D_POSIX_C_SOURCE=200112L -D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1
LIBS=-lblas

# primary target
cc4s: bin/Cc4s

OBJECTS=obj/Options.o obj/Chi.o obj/CoulombIntegrals.o obj/Amplitudes.o obj/Cc4s.o
# dependencies
obj/Cc4s.o: obj/Options.o obj/Chi.o obj/CoulombIntegrals.o obj/Amplitudes.o

# create directories of not present

clean:
	rm -rf bin/*
	rm -rf obj/*

# compile object files
obj/%.o: src/%.cxx src/%.hpp
	${CXX} ${COPTIONS} ${OPTIMIZE} -c src/$*.cxx -o $@ -I${CTF}/include

# compile and link executable
bin/%: obj/%.o
	${CXX} ${COPTIONS} ${OPTIMIZE} ${OBJECTS} -o $@ -I${CTF}/include/ -L${CTF}/lib -lctf ${LIBS}


