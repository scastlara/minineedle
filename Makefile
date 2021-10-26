PYTHON=python3.9
VENV_TEST=venv/venv_test


$(VENV_TEST): requirements/testing.txt
	$(PYTHON) -m venv venv
	venv/bin/pip install -r requirements/testing.txt
	touch $@

venv: $(VENV_TEST)

lint: venv
	venv/bin/isort --diff --check minineedle tests
	venv/bin/black --diff --check minineedle tests
	venv/bin/mypy minineedle

unittests: venv
	venv/bin/pytest

format:
	venv/bin/isort minineedle tests
	venv/bin/black minineedle tests

test: lint unittests

