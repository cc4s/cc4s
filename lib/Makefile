CONFIG ?= gxx

$(info )
$(info Using configuration $(CONFIG))
$(info )

include config.${CONFIG}
include Libraries

.DEFAULT_GOAL := all

all: ${LIBRARIES}

clean:
	rm -rf build/$(CONFIG)/
