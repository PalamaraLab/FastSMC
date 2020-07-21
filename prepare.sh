#!/bin/bash

cd TOOLS/PREPARE_DECODING || exit

# demographics file
demoFile=../../FILES/CEU.demo

# Time discretization. (See FILES/DISC/*)
states=30-100-2000
discFile=../../FILES/DISC/${states}.disc

# number of CSFS individuals. Values >= 50 should be sufficient, 300 is good. Needs to be <= number of (haploid) samples in the data set
CSFSsamples=100

# root of file to be analyzed
outFile=../../FILES/EXAMPLE/exampleFile.n100.array

# frequency file
freqFile=../../FILES/UKBB.frq

# get csfs
echo "Building csfs"
python3 get_csfs.py -D ${demoFile} -d ${discFile} -n ${CSFSsamples} -o ${outFile} > /dev/null || exit 1
echo "Finished building csfs"

# precompute decoding quantities
echo "Precomputing decoding quantities"
java -jar ASMCprepareDecoding.jar \
	-D ${demoFile} \
	-d ${discFile} \
	-n ${CSFSsamples} \
	-F ${freqFile} \
	-C ${outFile}.csfs \
	-o ${outFile} || exit 1
echo "Finished precomputing decoding quantities"
