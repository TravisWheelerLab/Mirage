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
BUILD_DIR=build

HSI_DIR=${DEPS_DIR}/${HSI_NAME}
SPALN_DIR=${DEPS_DIR}/${SPALN_NAME}
BLAT_DIR=${DEPS_DIR}/${BLAT_NAME}
TBLASTN_DIR=${DEPS_DIR}/${TBLASTN_NAME}

HSI_TAR=${HSI_DIR}.tgz
SPALN_TAR=${SPALN_DIR}.tgz
BLAT_TAR=${BLAT_DIR}.tgz
TBLASTN_TAR=${TBLASTN_DIR}.tgz

if [ ! -d ${BUILD_DIR} ];
then
    echo "No build found, see README"
    exit 1
fi

# Add hsi to build
if [ ! -d ${BUILD_DIR}/${HSI_NAME} ];
then
    tar -xf $HSI_TAR -C $DEPS_DIR
    cd $HSI_DIR && cmake . && make && cd -
    mv $HSI_DIR ${BUILD_DIR}/${HSI_BASE}
else
    echo "Skipping hsi"
fi

# Add spaln to build
if [ ! -d ${BUILD_DIR}/${SPALN_NAME} ];
then
    tar -xf $SPALN_TAR -C ${DEPS_DIR}
    cd $SPALN_DIR/src && ./configure && make && cd -
    mv $SPALN_DIR ${BUILD_DIR}/${SPALN_BASE}
else
    echo "Skipping spaln"
fi

# Add blat to build
if [ ! -d ${BUILD_DIR}/${BLAT_NAME} ];
then
    tar -xf $BLAT_TAR -C ${DEPS_DIR}
    mv $BLAT_DIR ${BUILD_DIR}/${BLAT_BASE}
else
    echo "Skipping blat"
fi

# Add tblastn to build
if [ ! -d ${BUILD_DIR}/${TBLASTN_NAME} ];
then
    tar -xf $TBLASTN_TAR -C ${DEPS_DIR}
    mv $TBLASTN_DIR ${BUILD_DIR}/${TBLASTN_BASE}
else
    echo "Skipping tblastn"
fi
