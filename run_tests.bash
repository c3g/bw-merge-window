#!/usr/bin/env bash

poetry install
poetry run pytest -svv --cov=bw_merge_window --cov-branch
