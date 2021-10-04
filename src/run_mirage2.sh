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

if [ -e "$path_only/build/Mirage2.pl" ];
then
    MIRAGECMD="$path_only/build/Mirage2.pl $@"
    $MIRAGECMD
elif [ -e "$path_only/Mirage2.pl" ];
then
     MIRAGECMD="$path_only/Mirage2.pl $@"
     $MIRAGECMD
else
    echo "\n  Failed to locate Mirage2.pl\n\n"
fi
