#!/bin/bash

#echo "Downloading and building Boost"
#sh getBoost.sh

echo "Building ASMC"
cd ASMC_SRC/SRC && make clean && make && cd ../BIN && ln -s $(pwd)/ASMC ../../
