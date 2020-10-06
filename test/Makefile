PYTHON ?= python3
ENV ?= env
PYENV = $(ENV)/bin/python
PIPENV = $(ENV)/bin/pip
TESTIS = $(ENV)/bin/testis
CONFIG ?= gxx-debug
CC4S_PATH ?= $(shell readlink -f ../build/$(CONFIG)/bin/Cc4s)
export CC4S_PATH

test: $(TESTIS)
	. $(ENV)/bin/activate && $< *

$(ENV)/.piped: $(ENV) requirements.txt
	$(PIPENV) install -r requirements.txt
	touch $@

$(ENV):
	virtualenv -p $(PYTHON) env

$(TESTIS): $(ENV)/.piped
	$(PYENV) setup.py develop

.PHONY: test
