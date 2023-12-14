#!/usr/bin/env bash
mkdir -p htmlcov
docker run --platform linux/x86_64 -v ./tests/out:/test/tests/out -v ./htmlcov:/test/htmlcov bwmw-test
docker container prune -f
