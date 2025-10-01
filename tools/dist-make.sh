#!/usr/bin/env bash
# Script to create a distribution tarball

set -ex

YAML_CPP_GIT_TAG=$(grep -m1 YAML_CPP_GIT_TAG CMakeLists.txt | awk '{print $2}' | tr -d "\")")
CTF_GIT_TAG=$(grep -m1 CTF_GIT_TAG CMakeLists.txt | awk '{print $2}' | tr -d "\")")

echo ${YAML_CPP_GIT_TAG}
echo ${CTF_GIT_TAG}
# Define the version
VERSION=$(git describe --all --dirty --long | tr "[/]" "_")
PACKAGE_NAME="cc4s-${VERSION}"

# Create the distribution directory
DIST_ROOT=dist
DIST_DIR="${DIST_ROOT}/${PACKAGE_NAME}"
mkdir -p "$DIST_DIR"

# Remove any existing tarball
rm -f "$DIST_DIR/$PACKAGE_NAME.tar.gz"

# Copy the source code
git archive --format=tar HEAD | tar -x -C "${DIST_DIR}"

##################################################
# Third party
#
mkdir -p ${DIST_DIR}/thirdparty

CTF_DIR=${DIST_DIR}/thirdparty/ctf
git clone https://gitlab.cc4s.org/cc4s/ctf ${CTF_DIR}
(
  cd ${CTF_DIR}
  git checkout ${CTF_GIT_TAG}
)

YAML_CPP_DIR=${DIST_DIR}/thirdparty/yaml-cpp
git clone https://gitlab.cc4s.org/cc4s/yaml-cpp ${YAML_CPP_DIR}
(
  cd ${YAML_CPP_DIR}
  git checkout ${YAML_CPP_GIT_TAG}
)

PYYAML_DIR=${DIST_DIR}/thirdparty/pyyaml
git clone https://github.com/yaml/pyyaml ${PYYAML_DIR}

##################################################
# testis-store
#
(
  cd ${DIST_DIR}/test
  make data
)

# Create the tarball
tar -czvf "${DIST_ROOT}/$PACKAGE_NAME.tar.gz" "${DIST_DIR}"

echo "Distribution tarball created: $DIST_DIR/$PACKAGE_NAME.tar.gz"
