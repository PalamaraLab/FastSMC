# script to build CSFS
# input: demographicModelFile discretizationFile numberOfTotalSamples outputRoot

demographicModelFile=$1
discretizationFile=$2
numberOfTotalSamples=$3
outputRoot=$4

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python3 $DIR/get_csfs.py -D $demographicModelFile -d $discretizationFile -n $numberOfTotalSamples -o $outputRoot > /dev/null
