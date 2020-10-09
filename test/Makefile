PYTHON ?= python3
TESTIS = ./testis.py
CONFIG ?= gxx-debug

include config.$(CONFIG)


ifndef CC4S_PATH
$(error DEAR USER: Please define the CC4S_PATH environment variable)
else
export CC4S_PATH
endif

$(info CC4S_PATH = $(CC4S_PATH))

test:
	$(PYTHON) $(TESTIS) -r . -n $(OUTPUT_NAME)

check:
	$(PYTHON) $(TESTIS) -c . -n $(OUTPUT_NAME)

all:
	$(PYTHON) $(TESTIS) -rc . -n $(OUTPUT_NAME)

.PHONY: test check all
