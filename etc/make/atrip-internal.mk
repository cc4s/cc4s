ATRIP_GIT_REPOSITORY ?= https://github.com/alejandrogallo/atrip
ATRIP_SRC_FILES = atrip/$(ATRIP_COMMIT)/src/atrip/Atrip.cxx

$(ATRIP_SRC_PATH):
	mkdir -p $@
	git clone $(ATRIP_GIT_REPOSITORY) $@
	cd $@ && git checkout $(ATRIP_COMMIT)

INTERNAL_DEPENDENCIES += $(ATRIP_SRC_PATH)
FETCH_DEPENDENCIES += $(ATRIP_SRC_PATH)
REQUIRED_MAKEVARS += ATRIP_COMMIT ATRIP_GIT_REPOSITORY ATRIP_SRC_PATH
$(ATRIP_SRC_FILES): $(ATRIP_SRC_PATH)
SRC_FILES += $(ATRIP_SRC_FILES)
INCLUDE_FLAGS += -I$(ATRIP_SRC_PATH)/include
