CONFIG ?= gcc
-include config.mk
CONFIG_PATH ?= etc/config/${CONFIG}.mk
include ${CONFIG_PATH}
include Sources.mk

# Check that reqiured makefile variables are indeed defined
define assert-vardefined
$(if $(strip $($(1))),,$(error $(1) must be defined))
endef
$(foreach _a,$(REQUIRED_MAKEVARS),$(call assert-vardefined,$(_a)))

# goals:
.DEFAULT_GOAL := cc4s
.PHONY: test cc4s clean extern

BUILD_PATH        = build/$(CONFIG)
OBJ_PATH          = $(BUILD_PATH)/obj
BIN_PATH          = $(BUILD_PATH)/bin
OBJ_FILES         = $(patsubst %.cxx,$(OBJ_PATH)/%.o,$(SRC_FILES))
DEP_FILES         = $(patsubst %.cxx,$(OBJ_PATH)/%.d,$(SRC_FILES))
TESTS_OBJECTS     = $(patsubst %.cxx,$(OBJ_PATH)/%.o,$(TEST_SRC_FILES))
# retrieve build environment
VERSION          := $(shell git describe --all --dirty --long)
DATE             := $(shell git log -1 --format="%cd")
COMPILER_VERSION := $(shell ${CXX} --version | head -n 1)


# This is a trick just to make sure that the dependencies are built
EXTERN_DONE_FILE = $(BUILD_PATH)/extern-built
$(EXTERN_DONE_FILE): $(EXTERNAL_DEPENDENCIES)
	@mkdir -p $(@D)
	@touch $@
extern: $(EXTERN_DONE_FILE)


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),clean-all)
ifneq ($(MAKECMDGOALS),extern)
ifneq ($(CREATE_EXTERN),FALSE)
ifneq ($(wildcard $(EXTERN_DONE_FILE)), $(EXTERN_DONE_FILE))
$(info )
$(info )
$(info You have not built the dependencies for $(CONFIG), please build them)
$(info with:     $(MAKE) extern CONFIG=$(CONFIG))
$(info )
$(info )
$(error exiting now...)
endif
# try to include dependency files in order to trigger their creation
-include $(DEP_FILES)
endif
endif
endif
endif

cc4s: $(BIN_PATH)/${CC4S_TARGET}

clean:
	rm -rf $(OBJ_PATH)
	rm -rf $(BIN_PATH)

clean-all: clean
	@$(MAKE) CREATE_EXTERN=FALSE $(patsubst %,%-clean,$(EXTERNAL_DEPENDENCIES))

test:
	$(MAKE) -C $@

unit-test: $(BIN_PATH)/Test

# generate documentation
doc:
	doxygen

# copy binary to installation directory
install: $(BIN_PATH)/$(CC4S_TARGET)
	mkdir -p ${CC4S_INSTALL}
	cp $< ${CC4S_INSTALL}

# build dependencies only
depend: $(DEP_FILES)

# keep intermediate files
.PRECIOUS: ${OBJ_FILES} ${DEP_FILES}

# create a dependency for object file
$(OBJ_PATH)/%.d: src/%.cxx
	$(info [DEP] $@)
	mkdir -p $(dir $@)
	${CXX} -MM ${CXXFLAGS} ${INCLUDE_FLAGS} -c src/$*.cxx | \
	  sed 's#[^ :]*\.o[ :]*#$(OBJ_PATH)/$*.o $@: #g' > $@

# compile an object file
$(OBJ_PATH)/%.o: $(OBJ_PATH)/%.d
	$(info [OBJ] $@)
	mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${INCLUDE_FLAGS} -c src/$*.cxx -o $@

# compile and link executable
$(BIN_PATH)/${CC4S_TARGET}: ${OBJ_FILES}
	$(info [BIN] $@)
	mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${OBJ_FILES} ${LDFLAGS} -o $@

# compile and link test executable
$(BIN_PATH)/Test: ${OBJ_FILES} $(TESTS_OBJECTS)
	mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${OBJ_FILES} $(TESTS_OBJECTS) ${LDFLAGS} -o $@
