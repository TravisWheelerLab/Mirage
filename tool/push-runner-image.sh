#!/usr/bin/env sh

set -e

# TODO: Update this to read the actual version once we have one
MIRAGE_VERSION=2.0.0-alpha

docker push traviswheelerlab/mirage:${MIRAGE_VERSION}
docker push traviswheelerlab/mirage:latest
