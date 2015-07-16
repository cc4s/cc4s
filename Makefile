# location of the Cyclops Tensor Framework library
CTF=../ctf
#TODO: use configuration files
CXX=mpicxx
OPTIMIZE=-O3
COPTIONS=-std=c++0x -fopenmp -Wall -D_POSIX_C_SOURCE=200112L -D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1
LIBS=-lblas

# primary target
cc4s: bin/cc4s

OBJECTS=obj/CoulombIntegrals.o obj/Chi.o obj/cc4s.o
# dependencies
obj/cc4s.o: obj/CoulombIntegrals.o obj/Chi.o

# create directories of not present

clean:
	rm -rf bin/*
	rm -rf obj/*

# compile object files
obj/%.o: src/%.cxx
	${CXX} ${COPTIONS} ${OPTIMIZE} -c src/$*.cxx -o $@ -I${CTF}/include

# compile and link executable
bin/%: obj/%.o
	${CXX} ${COPTIONS} ${OPTIMIZE} ${OBJECTS} -o $@ -I${CTF}/include/ -L${CTF}/lib -lctf ${LIBS}


