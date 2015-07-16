# location of the Cyclops Tensor Framework library
CTF=../ctf
#TODO: use configuration files
CXX=mpicxx
OPTIMIZE=-O3
COPTIONS=-std=c++0x -fopenmp -Wall -D_POSIX_C_SOURCE=200112L -D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1
LIBS=-lblas

# primary target
cc4s: bin/cc4s

# dependencies
bin/cc4s:

# create directories of not present
obj:
	mkdir -p obj
bin:
	mkdir -p bin

# compile object files
obj/%: obj src/%.cxx
	${CXX} ${COPTIONS} ${OPTIMIZE} src/$*.cxx -o $@ -I${CTF}/include

# compile and link executable
bin/%: bin src/%.cxx $(CTF)/lib/libctf.a 
	${CXX} ${COPTIONS} ${OPTMIZE} src/$*.cxx -o $@ -I${CTF}/include/ -L${CTF}/lib -lctf ${LIBS}


