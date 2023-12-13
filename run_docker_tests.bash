#!/usr/bin/env bash
docker run --platform linux/x86_64 -v ./tests/out:/test/tests/out bwmw-test
docker container prune -f
