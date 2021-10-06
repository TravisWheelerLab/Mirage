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

if [ -e "$path_only/build/Diviner.pl" ];
then
    DIVINERCMD="$path_only/build/Diviner.pl $@"
    $DIVINERCMD
elif [ -e "$path_only/Diviner.pl" ];
then
     DIVINERCMD="$path_only/Diviner.pl $@"
     $DIVINERCMD
else
    echo "\n  Failed to locate Diviner.pl\n\n"
fi

