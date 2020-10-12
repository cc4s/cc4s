YAML_SRC_PATH = $(abspath lib/src/yaml-cpp/$(YAML_COMMIT))
YAML_STATIC_LIB = $(YAML_BUILD_PATH)/libyaml-cpp.a
YAML_CMAKE = $(YAML_SRC_PATH)/CMakeLists.txt
CMAKE ?= cmake

$(YAML_CMAKE):
	mkdir -p $(@D)
	git clone https://github.com/jbeder/yaml-cpp.git $(@D)
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

IN_PROJECT_DEPENDENCIES += yaml
