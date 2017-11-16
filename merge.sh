# first input argument: root of file to be analyzed
file=$1
# number of batches to be merged
N=$2

# merge output of all jobs
java -jar MergePosteriorSums.jar \
        --fileRoot $file \
        --numJobs $N \
        --out $file
# remove individual batches
rm $file.*-*.??.sumOverPairs.gz
