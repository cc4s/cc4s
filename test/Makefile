PYTHON ?= python3
ENV_FOLDER ?= env
PYENV = $(ENV_FOLDER)/bin/python
PIPENV = $(ENV_FOLDER)/bin/pip
TESTIS = $(ENV_FOLDER)/bin/testis
CONFIG ?= gxx-debug
CC4S_PATH ?= $(shell readlink -f ../build/$(CONFIG)/bin/Cc4s)
export CC4S_PATH

test: $(TESTIS)
	. $(ENV_FOLDER)/bin/activate && $< *

$(ENV_FOLDER)/.piped: $(ENV_FOLDER) requirements.txt
	$(PIPENV) install -r requirements.txt
	touch $@

$(ENV_FOLDER):
	virtualenv -p $(PYTHON) env

$(TESTIS): $(ENV_FOLDER)/.piped
	$(PYENV) setup.py develop

.PHONY: test
