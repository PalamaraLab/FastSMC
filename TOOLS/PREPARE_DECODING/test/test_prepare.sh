#!/bin/bash

cd .. || exit 1

# demographics file
demoFile=../../FILES/CEU.demo

# Time discretization. (See FILES/DISC/*)
states=30-100-2000
discFile=../../FILES/DISC/${states}.disc

# number of CSFS individuals. Values >= 50 should be sufficient, 300 is good. Needs to be <= number of (haploid) samples in the data set
CSFSsamples=50

# root of file to be analyzed
outFile=test/test_output

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

python3 test/test_file_contents.py test/test_output.csfs test/reference_output.csfs || exit 1
python3 test/test_file_contents.py test/test_output.intervalsInfo test/reference_output.intervalsInfo || exit 1
python3 test/test_file_contents.py test/test_output.decodingQuantities.gz test/reference_output.decodingQuantities.gz || exit 1

echo "Successfully checked files"

rm test/test_output.csfs
rm test/test_output.intervalsInfo
rm test/test_output.decodingQuantities.gz
