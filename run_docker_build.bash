#!/usr/bin/env bash
docker build --platform linux/x86_64 . -f test.Dockerfile -t bwmw-test
