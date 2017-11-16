# this is needed on my machine to load python3
module load gcc/6.2.0

# first input argument: root of file to be analyzed
file=$1
# change this to modify the number of hidden coalescent states. They're precomputed in DISC for convenience
states=30-100-2000
# number of CSFS individuals. Values >= 50 should be sufficient. Needs to be <= number of (haploid) samples in the data set
CSFSsamples=100

### prepare decoding (done once, only depends on demographic model and number of states)
# first get csfs
echo "Building csfs"
python3 get_csfs.py -D FILES/CEU.demo -d FILES/DISC/$states.disc -n $CSFSsamples -o $file > /dev/null
# then correct csfs for array ascertainment, precompute quantities used for linear-time decoding
java -jar PrepareCoalescentDecoding.jar -D FILES/CEU.demo -d FILES/DISC/$states.disc -n $CSFSsamples -F FILES/chrAll.frq -C $file.csfs -o $file
