#!/bin/bash

# homebrew packages
brew update > /dev/null
brew install libomp

# install fastjet
git clone https://gitlab.com/pkomiske/fastjet.git
cd fastjet
git submodule update --init --recursive
autoreconf -i
export CXXFLAGS=-std=c++14
./configure --prefix=/usr/local --enable-pyext --enable-cgal --enable-cgal-header-only --disable-monolithic --disable-allplugins --disable-debug PYTHON=python3 PYTHON_CONFIG=python3-config
make -j2 install
cd ..
