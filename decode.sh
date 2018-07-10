#!/bin/bash

# first input argument: root of file to be analyzed
file=FILES/EXAMPLE/exampleFile.n100.array
# second input argument: ASMC mode (either "array" or "sequence")
mode=array
# number of batches to be analyzed. This can be parallelized
N=5

if [ "$mode" != "array" ] && [ "$mode" != "sequence" ]; then
        echo "ERROR: unrecognized mode: "$mode" please use one of {array, sequence}" > /dev/stderr;
        exit;
fi

### run analysis. This can be parallelized, here just using a loop
for i in $(seq 1 $N); do
	./ASMC --hapsFileRoot $file \
		--decodingQuantFile $file.decodingQuantities.gz \
	    --mode $mode \
		--majorMinorPosteriorSums \
		--jobs $N \
		--jobInd $i
done
