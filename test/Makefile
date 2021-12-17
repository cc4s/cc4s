-include config.mk
PYTHON ?= python3
CONFIG ?= gcc-debug
CC4S_PATH ?= $(realpath ../build/$(CONFIG)/bin/Cc4s)
CC4S_RUN ?= $(CC4S_PATH)
ifneq ($(MAKECMDGOALS),data)
ifneq ($(MAKECMDGOALS),pyyaml)
OUTPUT_NAME ?= $(CONFIG)-$(shell md5sum "$(CC4S_PATH)" | cut -d ' ' -f1).testis
endif
endif
TAGS ?= essential
TESTMATCH ?= .*

TESTIS = ./testis.py
TEST_FOLDERS = $(shell find . -name test.json | xargs dirname)
TESTIS_LINKS = $(patsubst %,%/$(TESTIS),$(TEST_FOLDERS))
BINARY_LINKS = $(patsubst %,%/$(OUTPUT_NAME)/Cc4s,$(TEST_FOLDERS))

-include config.$(CONFIG)

link: $(TESTIS_LINKS) $(BINARY_LINKS)
$(TESTIS_LINKS): $(TESTIS)
	ln -rs $< $@

$(BINARY_LINKS):
	mkdir -p $(@D)
	ln -rs $(CC4S_PATH) $@

ifneq ($(MAKECMDGOALS),data)
ifneq ($(MAKECMDGOALS),pyyaml)
ifeq ($(CC4S_RUN),)
$(error DEAR USER: Please define the CC4S_RUN environment variable)
else
export CC4S_PATH CC4S_RUN
endif
endif
endif

$(info CC4S_PATH = $(CC4S_PATH))
$(info CC4S_RUN = $(CC4S_RUN))

test: link
test:
	$(PYTHON) $(TESTIS) -rc . -n $(OUTPUT_NAME) \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

check: link
check:
	$(PYTHON) $(TESTIS) -c . -n $(OUTPUT_NAME) \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

run: link
run:
	$(PYTHON) $(TESTIS) -r . -n $(OUTPUT_NAME) \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

data:
	$(PYTHON) $(TESTIS) . \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

clean:
	-@find . -name '*.testis'   -type d -exec rm -vr {} \;
	-@find . -name testis-store -type d -exec rm -vr {} \;
	-@find . -name __pycache__  -type d -exec rm -vr {} \;
	-@rm -v $(TESTIS_LINKS)

clean-all: clean
	rm -vrf lib
	rm -vfr pyyaml

gdb: CC4S_RUN="gdb $(CC4S_PATH)"

pyyaml: lib/python/yaml
lib/python/yaml:
	git clone https://github.com/yaml/pyyaml
	mkdir -p $(@D)
	ln -rs pyyaml/lib/yaml/ $@
	

.PHONY: test check all gdb clean data pyyaml clean-all
