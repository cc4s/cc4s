TESTIS = ./testis.py
PYTHON ?= python3
CONFIG ?= gcc-debug
CC4S_PATH ?= $(realpath ../build/$(CONFIG)/bin/Cc4s)
CC4S_RUN ?= $(CC4S_PATH)
OUTPUT_NAME ?= $(CONFIG)-$(shell md5sum $(CC4S_PATH) | cut -d ' ' -f1)
TAGS =
TESTMATCH ?= .*

-include config.$(CONFIG)


ifeq ($(CC4S_RUN),)
$(error DEAR USER: Please define the CC4S_RUN environment variable)
else
export CC4S_PATH CC4S_RUN
endif

$(info CC4S_PATH = $(CC4S_PATH))
$(info CC4S_RUN = $(CC4S_RUN))

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
