#!/usr/bin/env sh

set -e

# TODO: Update this to read the actual version once we have one
MIRAGE_VERSION=2.0.0-alpha

docker build -f Dockerfile_run \
    -t traviswheelerlab/mirage:${MIRAGE_VERSION} \
    -t traviswheelerlab/mirage:latest \
    $@ \
    .
