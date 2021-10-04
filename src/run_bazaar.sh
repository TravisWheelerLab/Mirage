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

if [ -e "$path_only/build/Bazaar.pl" ];
then
    BAZAARCMD="$path_only/build/Bazaar.pl $@"
    $BAZAARCMD
elif [ -e "$path_only/Bazaar.pl" ];
then
     BAZAARCMD="$path_only/Bazaar.pl $@"
     $BAZAARCMD
else
    echo "\n  Failed to locate Bazaar.pl\n\n"
fi

