#!/bin/bash

cd ASMC_SRC/BOOST

echo "Downloading Boost"
wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz

echo "Extracting Boost (will take a while)..."
tar -xf boost_1_67_0.tar.gz

echo "Building required Boost libraries"
cd boost_1_67_0
./bootstrap.sh --with-libraries=program_options,iostreams
./b2
