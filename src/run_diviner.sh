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

if [ -e "$path_only/build/Oasis.pl" ];
then
    OASISCMD="$path_only/build/Oasis.pl $@"
    $OASISCMD
elif [ -e "$path_only/Oasis.pl" ];
then
     OASISCMD="$path_only/Oasis.pl $@"
     $OASISCMD
else
    echo "\n  Failed to locate Oasis.pl\n\n"
fi

