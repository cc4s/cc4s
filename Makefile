# default configuration, override with
#   make all CONFIG=(icc|gxx|icc-debug|icc-local|icc-debug-local)
CONFIG ?= gxx

include etc/config/${CONFIG}.mk
include Objects.mk

ifneq ($(IS_CLEANING),TRUE)
	# include created dependencies
	# if the makefile is compiling
	-include ${OBJECTS:.o=.d}
	-include build/${CONFIG}/obj/main/${CC4S_TARGET}.d
endif

# goals:
.DEFAULT_GOAL := all

deps: IS_CLEANING=TRUE
deps: $(IN_PROJECT_DEPENDENCIES)

clean: IS_CLEANING=TRUE
clean:
	rm -rf build/$(CONFIG)

# primary target
all: cc4s
cc4s: build/${CONFIG}/bin/${CC4S_TARGET}

.PHONY: test wiki cc4s all clean deps
test:
	$(MAKE) -C $@

unit-test: build/${CONFIG}/bin/Test

# generate documentation
doc:
	doxygen

wiki:
	bash utils/extract.sh -R -d wiki/dist -b wiki/build -p src test.wiki

# copy binary to installation directory
install: build/${CONFIG}/bin/${CC4S_TARGET}
	mkdir -p ${CC4S_INSTALL}
	cp build/${CONFIG}/bin/${CC4S_TARGET} ${CC4S_INSTALL}

# build dependencies only
depend: ${OBJECTS:.o=.d} build/${CONFIG}/obj/main/${CC4S_TARGET}.d

# retrieve build environment
VERSION:=$(shell git describe --all --dirty --long)
DATE:=$(shell git log -1 --format="%cd")
COMPILER_VERSION:=$(shell ${CXX} --version | head -n 1)

# add build environment specifics to INCLUDE_FLAGS and to CC4S_OPTIONS
INCLUDE_FLAGS += -Isrc/main
CC4S_OPTIONS += -D_POSIX_C_SOURCE=200112L \
-D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1 -DCC4S_VERSION=\"${VERSION}\" \
"-DCC4S_DATE=\"${DATE}\"" \
"-DCOMPILER_VERSION=\"${COMPILER_VERSION}\""


# create a dependency for object file
build/${CONFIG}/obj/%.d: src/%.cxx
	mkdir -p $(dir $@)
	${CXX} -MM ${CC4S_OPTIONS} ${INCLUDE_FLAGS} -c src/$*.cxx | \
	  sed 's#[^ :]*\.o[ :]*#build/${CONFIG}/obj/$*.o $@: #g' > $@

# keep dependency files
.PRECIOUS: build/${CONFIG}/obj/%.o ${OBJECTS}

# compile an object file
build/${CONFIG}/obj/%.o: build/${CONFIG}/obj/%.d
	mkdir -p $(dir $@)
	${CXX} ${CC4S_OPTIONS} ${OPTIMIZE} ${INCLUDE_FLAGS} -c src/$*.cxx -o $@

# keep object files
.PRECIOUS: build/${CONFIG}/obj/%.o ${OBJECTS}

# compile and link executable
build/${CONFIG}/bin/%: build/${CONFIG}/obj/main/%.o ${OBJECTS}
	mkdir -p $(dir $@)
	${CXX} ${CC4S_OPTIONS} ${OPTIMIZE} ${OBJECTS} build/${CONFIG}/obj/main/${CC4S_TARGET}.o ${INCLUDE_FLAGS} ${LINK_LIBS} -o $@

# compile and link test executable
build/${CONFIG}/bin/Test: ${OBJECTS} $(TESTS_OBJECTS)
	mkdir -p $(dir $@)
	${CXX} ${CC4S_OPTIONS} ${OPTIMIZE} ${OBJECTS} $(TESTS_OBJECTS) ${INCLUDE_FLAGS} ${LINK_LIBS} -o $@
