#!/bin/bash

realpath_ours()
{
    OURPWD=$PWD
    cd "$(dirname "$1")"
    LINK=$(readlink "$basename "$1")")
    while [ "$LINK" ]; do
	cd "$(dirname "$LINK")"
	LINK=$(readlink "$(basename "$1")")
    done
    REALPATH="$PWD/$(basename "$1")"
    cd "$OURPWD"
    echo "$REALPATH"
}

path_only=$(dirname $(realpath_ours "$0"))
MIRAGECMD="$path_only/src/Mirage2.pl $@"
$MIRAGECMD