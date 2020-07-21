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

**This repository contains the code and example files for the ASMC program. A user manual can be found [here](http://www.palamaralab.org/software/ASMC), data and annotations from the paper can be found [here](http://www.palamaralab.org/data/ASMC).**

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

**C++ library only**
```bash
git clone https://github.com/OxfordRSE/ASMC.git

mkdir ASMC_BUILD_DIR && cd ASMC_BUILD_DIR
cmake /path/to/ASMC
cmake --build .
```

**C++ library and Python bindings**
```bash
git clone --recurse-submodules https://github.com/OxfordRSE/ASMC.git

cd ASMC
pip install .
```

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

## License

ASMC is distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on ASMC, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

If you use this software, please cite:

- P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. *Nature Genetics*, 2018.

