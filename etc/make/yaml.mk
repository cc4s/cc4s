YAML_STATIC_LIB = $(YAML_BUILD_PATH)/libyaml-cpp.a
YAML_CMAKE = $(YAML_SRC_PATH)/CMakeLists.txt
YAML_GIT_REPOSITORY ?= https://github.com/jbeder/yaml-cpp.git
CMAKE ?= cmake

$(YAML_CMAKE):
	mkdir -p $(@D)
	git clone  ${YAML_GIT_REPOSITORY} $(@D)
	cd $(@D) && git checkout $(YAML_COMMIT)

$(YAML_BUILD_PATH)/Makefile: $(YAML_CMAKE)
	mkdir -p $(YAML_BUILD_PATH)
	cd $(YAML_BUILD_PATH) && \
	$(CMAKE) -DCMAKE_INSTALL_PREFIX=$(YAML_BUILD_PATH) $(YAML_SRC_PATH)

$(YAML_STATIC_LIB): $(YAML_BUILD_PATH)/Makefile
	$(info Compiling $@)
	cd $(@D) && $(MAKE) install

.PHONY: yaml yaml-clean
yaml: $(YAML_STATIC_LIB)

yaml-clean:
	rm -rf $(YAML_BUILD_PATH)

EXTERNAL_DEPENDENCIES += yaml
INCLUDE_FLAGS += $(YAML_INCLUDE)
STATIC_LIBS += $(YAML_LDFLAGS)
REQUIRED_MAKEVARS += YAML_COMMIT YAML_BUILD_PATH YAML_GIT_REPOSITORY YAML_SRC_PATH
