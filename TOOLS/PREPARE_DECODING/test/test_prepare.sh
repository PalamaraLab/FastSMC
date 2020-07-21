#!/bin/bash

cd .. || exit 1

# demographics file
demoFile=../../FILES/CEU.demo

# Time discretization. (See FILES/DISC/*)
states=30-100-2000
discFile=../../FILES/DISC/${states}.disc

# number of CSFS individuals. Values >= 50 should be sufficient, 300 is good. Needs to be <= number of (haploid) samples in the data set
CSFSsamples=300

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

python3 test/test_file_contents.py test/test_output.csfs ../../FILES/DECODING_QUANTITIES/30-100-2000.csfs || exit 1
python3 test/test_file_contents.py test/test_output.intervalsInfo ../../FILES/DECODING_QUANTITIES/30-100-2000.intervalsInfo || exit 1
python3 test/test_file_contents.py test/test_output.decodingQuantities.gz ../../FILES/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz || exit 1

echo "Successfully checked files"

echo "another line" >> test/test_output.csfs
python3 test/test_file_contents.py test/test_output.csfs ../../FILES/DECODING_QUANTITIES/30-100-2000.csfs || exit 1

#rm test/test_output.csfs
#rm test/test_output.intervalsInfo
#rm test/test_output.decodingQuantities.gz
