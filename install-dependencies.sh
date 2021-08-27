#!/usr/bin/env sh

set -e

# Note: This script is run by the CMake build and there should be no reason to
# run it manually except for development and testing purposes.

# Note: when upgrading dependencies, be sure to update the base file names in
# the block below, where appropriate.

HSI_NAME=hsi-1.0.0
SPALN_NAME=spaln2.3.3
BLAT_NAME=blat
TBLASTN_NAME=tblastn

HSI_BASE=hsi
SPALN_BASE=spaln
BLAT_BASE=blat
TBLASTN_BASE=tblastn

DEPS_DIR=dependencies

HSI_DIR=${DEPS_DIR}/${HSI_NAME}
SPALN_DIR=${DEPS_DIR}/${SPALN_NAME}
BLAT_DIR=${DEPS_DIR}/${BLAT_NAME}
TBLASTN_DIR=${DEPS_DIR}/${TBLASTN_NAME}

HSI_TAR=${HSI_DIR}.tgz
SPALN_TAR=${SPALN_DIR}.tgz
BLAT_TAR=${BLAT_DIR}.tgz
TBLASTN_TAR=${TBLASTN_DIR}.tgz

if [ ! -d build/ ];
then
    echo "No build found, see README"
    exit 1
fi

# Add hsi to build
if [ ! -d build/${HSI_NAME} ];
then
    tar -xf $HSI_TAR -C $DEPS_DIR
    pushd $HSI_DIR
    cmake . && make
    popd
    mv $HSI_DIR build/${HSI_BASE}
else
    echo "Skipping hsi"
fi

# Add spaln to build
if [ ! -d build/${SPALN_NAME} ];
then
    tar -xf $SPALN_TAR -C dependencies/
    pushd $SPALN_DIR
    ./configure && make
    popd
    mv $SPALN_DIR build/${SPALN_BASE}
else
    echo "Skipping spaln"
fi

# Add blat to build
if [ ! -d build/${BLAT_NAME} ];
then
    tar -xf $BLAT_TAR -C dependencies/
    mv $BLAT_DIR build/${BLAT_BASE}
else
    echo "Skipping blat"
fi

# Add tblastn to build
if [ ! -d build/${TBLASTN_NAME} ];
then
    tar -xf $TBLASTN_TAR -C dependencies
    mv $TBLASTN_DIR build/${TBLASTN_BASE}
else
    echo "Skipping tblastn"
fi
