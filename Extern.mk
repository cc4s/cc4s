ifeq ($(HF), yes)
include etc/make/eigen.mk
include etc/make/libint.mk
endif

# base dependencies
include etc/make/yaml.mk
include etc/make/ctf.mk

ifeq ($(ATRIP),yes)
include etc/make/atrip-internal.mk
endif
