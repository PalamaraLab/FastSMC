[![Ubuntu unit tests](https://github.com/OxfordRSE/ASMC/workflows/Ubuntu%20unit/badge.svg)](https://github.com/OxfordRSE/ASMC/actions)
[![macOS unit tests](https://github.com/OxfordRSE/ASMC/workflows/macOS%20unit/badge.svg)](https://github.com/OxfordRSE/ASMC/actions)
[![Python tests](https://github.com/OxfordRSE/ASMC/workflows/Python%203.5%203.8/badge.svg)](https://github.com/OxfordRSE/ASMC/actions)
[![Regression test](https://github.com/OxfordRSE/ASMC/workflows/Regression%20test/badge.svg)](https://github.com/OxfordRSE/ASMC/actions)
[![Ubuntu asan](https://github.com/OxfordRSE/ASMC/workflows/Ubuntu%20asan/badge.svg)](https://github.com/OxfordRSE/ASMC/actions)
[![Ubuntu no sse/avx](https://github.com/OxfordRSE/ASMC/workflows/Ubuntu%20no%20sse/avx/badge.svg)](https://github.com/OxfordRSE/ASMC/actions)
[![codecov](https://codecov.io/gh/OxfordRSE/ASMC/branch/master/graph/badge.svg)](https://codecov.io/gh/OxfordRSE/ASMC)

```
    █████╗   ███████╗  ███╗   ███╗   ██████╗
   ██╔══██╗  ██╔════╝  ████╗ ████║  ██╔════╝
   ███████║  ███████╗  ██╔████╔██║  ██║     
   ██╔══██║  ╚════██║  ██║╚██╔╝██║  ██║     
   ██║  ██║  ███████║  ██║ ╚═╝ ██║  ╚██████╗
   ╚═╝  ╚═╝  ╚══════╝  ╚═╝     ╚═╝   ╚═════╝
```

The Ascertained Sequentially Markovian Coalescent is a method to efficiently estimate pairwise coalescence time along the genome. It can be run using SNP array or whole-genome sequencing (WGS) data.

**This repository contains code, installation instructions, and example files for the ASMC program. A user manual can be found [here](https://palamaralab.github.io/software/ASMC), data and annotations from the ASMC paper can be found [here](https://palamaralab.github.io/data/ASMC).**

## Installation

ASMC is regularly built and tested on Ubuntu and macOS.
It is a C++ library with optional Python bindings.

The ASMC C++ library requires:

 - A C++ compiler (C++14 or later)
 - CMake (3.12 or later)
 - Boost (1.62 or later)
 - Eigen (3.3.4 or later)

 The Python bindings additionally require:

 - Python (3.5 or later)
 - PyBind11 (distributed with ASMC as a submodule)

### Install dependencies

**Ubuntu (using the package manager)**
```bash
sudo apt install g++ cmake libboost-all-dev libeigen3-dev
```

**macOS (using homebrew and assuming Xcode is installed)**
```bash
brew install cmake boost eigen
```

### Getting and compiling ASMC

**C++ library only**:
```bash
git clone https://github.com/OxfordRSE/ASMC.git

mkdir ASMC/ASMC_BUILD_DIR && cd ASMC/ASMC_BUILD_DIR
cmake ..
cmake --build .
```

Note: you can locate the build directory outside the `ASMC` directory if you wish: just run `cmake /path/to/ASMC` from any directory you like.

**C++ library and Python bindings**
```bash
git clone --recurse-submodules https://github.com/OxfordRSE/ASMC.git

cd ASMC
pip install .
```

Note: the `--recurse-submodules` is important as PyBind11 is distributed as a submodule.
You will **not** get PyBind11 if you download the zip archive from GitHub.

### Decoding Quantities

To generate decoding quantities, several additional requirements are required.

**Ubuntu (using the package manager)**
```bash
sudo apt install libgmp-dev libmpfr-dev libgsl0-dev default-jdk jblas
```

**macOS (using homebrew and assuming cask is installed)**
```bash
brew install mpfr gmp gsl
brew cask install java
```

**Install python dependencies**
```bash
pip install cython numpy
pip install -r TOOLS/PREPARE_DECODING/requirements.txt
```

Basic functionality for generating decoding quantities can be seen in:

./prepare.sh


```
   ███████╗  █████╗  ███████╗ ████████╗ ███████╗ ███╗   ███╗  ██████╗
   ██╔════╝ ██╔══██╗ ██╔════╝ ╚══██╔══╝ ██╔════╝ ████╗ ████║ ██╔════╝
   █████╗   ███████║ ███████╗    ██║    ███████╗ ██╔████╔██║ ██║     
   ██╔══╝   ██╔══██║ ╚════██║    ██║    ╚════██║ ██║╚██╔╝██║ ██║     
   ██║      ██║  ██║ ███████║    ██║    ███████║ ██║ ╚═╝ ██║ ╚██████╗
   ╚═╝      ╚═╝  ╚═╝ ╚══════╝    ╚═╝    ╚══════╝ ╚═╝     ╚═╝  ╚═════╝
```

The Fast Sequentially Markovian Coalescent (FastSMC) algorithm is an extension to the ASMC algorithm, adding an identification step by hashing (currently using an improved version of the GERMLINE algorithm). FastSMC is an accurate method to detect Identical-By-Descent segments which enables estimating the time to most recent common ancestor for IBD individuals, and provides an estimate of uncertainty for detected IBD regions.

## Installation

FastSMC is compiled with ASMC, using the same instructions as above.

## Running FastSMC

You can run FastSMC as a C++ compiled executable or using Python (see below for examples).

**This document is not intended as an extensive guide, a more detailed user manual can be found [here](https://palamaralab.github.io/software/FastASMC), data and annotations from the FastSMC paper can be found [here](https://palamaralab.github.io/data/FastASMC).**

### Detailed command line options
See ASMC's documentation for parameters related to the validation step. Additional parameters related to the identification step are listed below. Note: default parameter values are likely to change in future versions.

```
  --inFileRoot                	Prefix of input files (.hap, .samples, .map).
                              	[mandatory]
  --decodingQuantFile         	Decoding quantities file.
                              	[mandatory]
  --outFileRoot               	Prefix of output file.
                              	[mandatory]
  --GERMLINE                  	Use of GERMLINE to pre-process IBD segments. If off, no identification step will be performed.
                              	[default 0/off]
  --min_m arg (=1)		Minimum match length (in cM).
				[default = 1.0]
  --time arg (=100)		Time threshold to define IBD in number of generations.
				[default = 100]
  --skip arg (=0)		Skip words with (seeds/samples) less than this value
				[default 0.0]
  --min_maf arg (=0)		Minimum minor allele frequency
				[default 0.0]
  --gap arg (=1)		Allowed gaps
                                [default 1]
  --max_seeds arg (=0)		Dynamic hash seed cutoff
				[default 0/off]
  --recall arg (=3)		Recall level from 0 to 3 (higher value means higher recall).
				[default = 3]
  --segmentLength		Output length in centimorgans of each IBD segment.
				[default 0/off]
  --perPairMAP			Output MAP age estimate for each IBD segment.
				[default 0/off]
  --perPairPosteriorMeans	Output posterior mean age estimate for each IBD segment.
				[default 0/off]
  --noConditionalAgeEstimates	Do not condition the age estimates on the TMRCA being between present time and t generations ago
				(where t is the time threshold).
				[default 0/off]
  --bin				Binary output
				[default off]
  --batchSize			Size of batches to be decoded.
				[default = 32]
```

Suggested optimal parameters for IBD detection within the past 25, 50, 100, 150 and 200 generations are provided in the FastSMC paper.

### Input file formats

Input files are provided to FastSMC with the --inFileRoot option. You may want to look at files in FILES/FASTSMC_EXAMPLE/* for examples of the file formats described below.

#### Phased haplotypes in Oxford haps/sample format (.hap/.hap.gz, .samples)
These files are provided in input to FastSMC. The file format explained [here](https://www.cog-genomics.org/plink/2.0/formats#haps). These files are output by phasing programs like Eagle and Shapeit.

#### Genetic map (.map)
The genetic map provided in input to FastSMC has 4 columns with format "Physical_position Recombination_rate Genetic_position Mutation_rate". Genetic positions are in centimorgans, physical positions are in bp. The map can be optionally compressed using gzip.

#### Decoding quantities (.decodingQuantities.gz)
See the instructions above to generate decoding quantities files and the ASMC manual [here](https://palamaralab.github.io/software/ASMC) for more details.

### Output format

FastSMC generates an .ibd.gz (or .bibd.gz if binary output) file in the specified location.
Each line corresponds to a pairwise shared segment, with the following fields:

	0. First individual's family identifier
	1. First individual identifier
	2. First individual haplotype identifier (1 or 2)
	3. Second individual's family identifier
	4. Second individual identifier
	5. Second individual haplotype identifier (1 or 2)
	6. Chromosome number
	7. Starting position of the IBD segment (inclusive)
	8. Ending position of the IBD segment (inclusive)
	9. (optional) Length in centimorgans of IBD segment
	10. IBD score
	11. (optional) Average mean posterior age estimate of the IBD segment
	12. (optional) Average MAP age estimate of the IBD segment

### Binary output

If you use the --bin option, FastSMC will generate a compressed binary (.bib.gz) output. This can be then converted to text format (see notebooks for an example).

### Examples using the Python bindings (Python)

FastSMC can also be run using Python. See the `notebooks` directory for an example.
There are two Jupyter notebooks:
 - a [minimal working example](notebooks/fastsmc-minimal.ipynb), where sensible defaults for parameters are chosen automatically
 - a [more detailed example](notebooks/fastsmc.ipynb) that demonstrates how to customise parameters, how to convert the binary file to text format, and how to analyse the output if it is too large to fit in memory.

### Example using the compiled FastSMC executable (C++)

Following the compilation instructions above will create an executable

```
ASMC_BUILD_DIR/FastSMC_exe
```

which can be used by providing command line arguments summarised below.

Either way of running FastSMC will run it on a simulated dataset with default parameters values.
An output file with IBD segments will be generated, and run time should be less than 4s.

## License

ASMC and FastSMC are distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on ASMC, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Reference

If you use this software, please cite the appropriate reference(s) below.

The ASMC algorithm and software were developed in
- P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. *Nature Genetics*, 2018.

The FastASMC algorithm and software were developed in
- J. Nait Saada, G. Kalantzis, D. Shyr, F. Cooper, M. Robinson, A. Gusev, P. F. Palamara. Identity-by-descent detection across 487,409 British samples reveals fine-scale population structure, evolutionary history, and trait associations. *Nature Communications*, in press.
