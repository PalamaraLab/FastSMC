# root of file to be analyzed
file=FILES/EXAMPLE/exampleFile.n100.array
# number of batches to be merged
N=5

# merge output of all jobs
java -jar TOOLS/MERGE_POSTERIORS/ASMCmergePosteriorSums.jar \
        --fileRoot $file \
        --numJobs $N \
	--norm \
        --out $file

# remove individual batches
rm $file.*-*.??.sumOverPairs.gz
