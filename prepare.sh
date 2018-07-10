#!/bin/bash

cd TOOLS/PREPARE_DECODING

# root of file to be analyzed
outFile=../../FILES/EXAMPLE/exampleFile.n100.array
# Time discretization. (See FILES/DISC/*)
states=30-100-2000
# number of CSFS individuals. Values >= 50 should be sufficient, 300 is good. Needs to be <= number of (haploid) samples in the data set
CSFSsamples=100

# get csfs
echo "Building csfs"
python3 get_csfs.py -D ../../FILES/CEU.demo -d ../../FILES/DISC/$states.disc -n $CSFSsamples -o $outFile > /dev/null
# precompute decoding quantities
java -jar ASMCprepareDecoding.jar \
	-D ../../FILES/CEU.demo \
	-d ../../FILES/DISC/$states.disc \
	-n $CSFSsamples \
	-F ../../FILES/UKBB.frq \
	-C $outFile.csfs \
	-o $outFile
