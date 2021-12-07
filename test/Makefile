-include config.mk
TESTIS = ./testis.py
PYTHON ?= python3
CONFIG ?= gcc-debug
CC4S_PATH ?= $(realpath ../build/$(CONFIG)/bin/Cc4s)
CC4S_RUN ?= $(CC4S_PATH)
OUTPUT_NAME ?= $(CONFIG)-$(shell md5sum $(CC4S_PATH) | cut -d ' ' -f1).testis
TAGS =
TESTMATCH ?= .*
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

ifeq ($(CC4S_RUN),)
$(error DEAR USER: Please define the CC4S_RUN environment variable)
else
export CC4S_PATH CC4S_RUN
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

clean:
	-@find . -name '*.testis'   -type d -exec rm -vr {} \;
	-@find . -name testis-store -type d -exec rm -vr {} \;
	-@find . -name __pycache__  -type d -exec rm -vr {} \;

gdb: CC4S_RUN="gdb $(CC4S_PATH)"

.PHONY: test check all gdb clean
