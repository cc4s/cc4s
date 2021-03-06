ATRIP_GIT_REPOSITORY ?= https://github.com/alejandrogallo/atrip
ATRIP_MAKEFILE = $(ATRIP_SRC_PATH)/Makefile
ATRIP_SOURCES_FILE = $(ATRIP_SRC_PATH)/Sources.mk
ATRIP_STATIC_LIB = $(ATRIP_BUILD_PATH)/lib/libatrip.a
ATRIP_CONFIG = lib

-include $(ATRIP_SOURCES_FILE)


$(ATRIP_SOURCES_FILE):
	mkdir -p $(@D)
	git clone $(ATRIP_GIT_REPOSITORY) $(@D)
	cd $(@D) && git checkout $(ATRIP_COMMIT)

$(ATRIP_STATIC_LIB): $(ATRIP_SOURCES_FILE)
	cd $(ATRIP_SRC_PATH) && \
		$(MAKE) CONFIG=$(ATRIP_CONFIG) \
						LDFLAGS="$(filter-out -latrip,$(LDFLAGS))" \
						CXXFLAGS="$(INCLUDE_FLAGS) -I$(ATRIP_SRC_PATH)/include -fPIC" \
						CXX=$(CXX) \
						CTF_STATIC_LIB:="$(CTF_STATIC_LIB)" \
						static
	cd $(ATRIP_SRC_PATH) && \
		$(MAKE) install CONFIG=$(ATRIP_CONFIG) PREFIX=$(ATRIP_BUILD_PATH)

.PHONY: atrip

$(ATRIP_STATIC_LIB): $(patsubst %,$(ATRIP_SRC_PATH)/%,$(ATRIP_SOURCES))
atrip: $(ATRIP_STATIC_LIB)

EXTERNAL_DEPENDENCIES += atrip
STATIC_LIBS += $(ATRIP_LDFLAGS)
INCLUDE_FLAGS += $(ATRIP_INCLUDE)
REQUIRED_MAKEVARS += ATRIP_COMMIT ATRIP_BUILD_PATH ATRIP_GIT_REPOSITORY ATRIP_SRC_PATH
