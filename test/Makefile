TESTIS = ./testis.py
PYTHON ?= python3
CONFIG ?= gcc-debug
CC4S_PATH ?= $(realpath ../build/$(CONFIG)/bin/Cc4s)
OUTPUT_NAME ?= $(CONFIG)-$(shell md5sum $(CC4S_PATH) | cut -d ' ' -f1)
TAGS =
TESTMATCH ?= .*

-include config.$(CONFIG)


ifeq ($(CC4S_PATH),)
$(error DEAR USER: Please define the CC4S_PATH environment variable)
else
export CC4S_PATH
endif

$(info CC4S_PATH = $(CC4S_PATH))

test:
	$(PYTHON) $(TESTIS) -rc . -n $(OUTPUT_NAME) \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

check:
	$(PYTHON) $(TESTIS) -c . -n $(OUTPUT_NAME) \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

run:
	$(PYTHON) $(TESTIS) -r . -n $(OUTPUT_NAME) \
		--rx "$(TESTMATCH)" \
		--tags "$(TAGS)"

.PHONY: test check all
