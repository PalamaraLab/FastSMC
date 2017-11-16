# first input argument: root of file to be analyzed
file=$1
# second input argument: ASMC mode (either "array" or "sequence")
mode=$2
# batch to be analyzed
n=$3
# number of batches to be analyzed. This can be parallelized
N=$4

if [ "$mode" != "array" ] && [ "$mode" != "sequence" ]; then
        echo "ERROR: unrecognized mode: "$mode" please use one of {array, sequence}" > /dev/stderr;
        exit;
fi

### run analysis. This can be parallelized, but is done sequentially here using a for loop
echo "Running batch"$n" of "$N
./ASMC --hapsFileRoot $file \
	--decodingQuantFile $file.decodingQuantities.gz \
        --mode $mode \
	--majorMinorPosteriorSums \
	--jobs $N \
	--jobInd $n
