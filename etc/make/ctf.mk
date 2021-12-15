CTF_CONFIG_FLAGS =
CTF_STATIC_LIB = $(CTF_BUILD_PATH)/lib/libctf.a
CTF_SHARED_LIB = $(CTF_BUILD_PATH)/lib/libctf.so
CTF_GIT_REPOSITORY ?= https://github.com/cyclops-community/ctf
CTF_CONFIGURE = $(CTF_SRC_PATH)/configure

$(CTF_CONFIGURE):
	mkdir -p $(@D)
	git clone $(CTF_GIT_REPOSITORY) $(@D)
	cd $(@D) && git checkout $(CTF_COMMIT)

$(CTF_BUILD_PATH)/Makefile: $(CTF_SRC_PATH)/configure
	mkdir -p $(CTF_BUILD_PATH)
	cd $(CTF_BUILD_PATH) && $(CTF_SRC_PATH)/configure $(CTF_CONFIG_FLAGS)

$(CTF_STATIC_LIB): $(CTF_BUILD_PATH)/Makefile
	$(info Compiling $@)
	cd $(CTF_BUILD_PATH) && $(MAKE)

.PHONY: ctf ctf-clean
ctf: $(CTF_STATIC_LIB)

ctf-clean:
	rm -rf $(CTF_BUILD_PATH)

EXTERNAL_DEPENDENCIES += ctf
FETCH_DEPENDENCIES += $(CTF_CONFIGURE)
STATIC_LIBS += $(CTF_LDFLAGS)
INCLUDE_FLAGS += $(CTF_INCLUDE)
REQUIRED_MAKEVARS += CTF_COMMIT CTF_BUILD_PATH CTF_GIT_REPOSITORY CTF_SRC_PATH
