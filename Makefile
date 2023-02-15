PYTHON?=python3.11
VENV_TEST=venv/venv_test
VENV=.venv

export PYTHONPATH=minineedle

$(VENV_TEST):
	POETRY_VIRTUALENVS_IN_PROJECT=1 poetry env use $(PYTHON)
	poetry install

venv: $(VENV_TEST)

lint: venv
	$(VENV)/bin/black --diff --check minineedle tests
	$(VENV)/bin/ruff check minineedle tests
	$(VENV)/bin/mypy minineedle tests

unit-tests: venv
	$(VENV)/bin/pytest

format:
	$(VENV)/bin/ruff minineedle tests --fix
	$(VENV)/bin/black minineedle tests

test: lint unit-tests
