```
    █████╗   ███████╗  ███╗   ███╗   ██████╗
   ██╔══██╗  ██╔════╝  ████╗ ████║  ██╔════╝
   ███████║  ███████╗  ██╔████╔██║  ██║     
   ██╔══██║  ╚════██║  ██║╚██╔╝██║  ██║     
   ██║  ██║  ███████║  ██║ ╚═╝ ██║  ╚██████╗
   ╚═╝  ╚═╝  ╚══════╝  ╚═╝     ╚═╝   ╚═════╝
```

The Ascertained Sequentially Markovian Coalescent is a method to efficiently estimate pairwise coalescence time along the genome. It can be run using SNP array or whole-genome sequencing (WGS) data. This page describes compiling and running the ASMC program.

### TL;DR

Download ASMC [here](https://github.com/pierpal/ASMC), edit `ASMC/SRC/Makefile` to point to your Boost library headers and binaries, or try to install them using `getBoost.sh`. Compile ASMC using `sh build.sh`. To compute pairwise coalescence times for the following files containing SNP array data for 150 phased diploid samples:

- `FILES/EXAMPLE/exampleFile.n100.array.hap.gz`
- `FILES/EXAMPLE/exampleFile.n100.array.samples`
- `FILES/EXAMPLE/exampleFile.n100.array.map.gz`

(see [here](#formats) for file formats), you can run the following ASMC command:
```
ASMC/ASMC \
        --decodingQuantFile FILES/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz \
        --hapsFileRoot FILES/EXAMPLE/exampleFile.n100.array \
        --posteriorSums \
        --outFileRoot FILES/EXAMPLE/exampleFile.n100.array
```

This will generate `FILES/EXAMPLE/exampleFile.n100.array.sumOverPairs.gz`, which contains a matrix of size SxD, where S is the number of sites in the data, and D is the number of discrete coalescence time intervals defined in `FILES/DISC/30-100-2000.disc`. The [*i*,*j*]-th entry of each matrix contains the sum of posterior coalescence probabilities for all samples at SNP *i* and time *j*. A more detailed example is described [here](#example).

### Table of Contents

- [Download, reference, license](#reference)
- [Dependencies](#dependencies)
- [Compiling](#compiling)
- [Running ASMC](#running)
- [Complete example](#example)
- [Detailed command line options](#detailed)
- [Tools and scripts](#tools)
- [Precomputed decoding quantities](#precomputed)
- [Change log](#log)

### Download, reference, license <a name="reference" />

ASMC can be downloaded [here](https://github.com/pierpal/ASMC). If you use this software, please cite:

- P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. *Nature Genetics*, 2018.

For any questions or comments on ASMC, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

ASMC is distributed under the GNU General Public License v3.0 (GPLv3).

### Dependencies <a name="dependencies" />

ASMC uses the following third party libraries.

- [SMC++](https://github.com/popgenmethods/smcpp) is needed if you would like to run `ASMCprepareDecoding` to analyze data with a non-European demographic model and/or your own time discretization or SNP allele frequencies ([this section](#precomputed) lists available precomputed models).
- [Boost](https://www.boost.org/) is needed to run the ASMC program. You can either download and install your own version or try using the `getBoost.sh` script to download and build Boost v1.67.0.

The following libraries are used by the `ASMCprepareDecoding` program (no need to install them, they are included in the `ASMCprepareDecoding.jar` file):

- [Argparse4j](https://argparse4j.github.io/) version 0.7.0.
- [Apache Commons Math](https://commons.apache.org/proper/commons-math/) version 3.3.
- [JBlas](http://jblas.org/) version 1.2.4.

### Compiling <a name="compiling" />

ASMC is composed of two main programs:

- `ASMCprepareDecoding`, written in Java, precomputes all the information needed to estimate pairwise coalescence time in SNP or WGS data.
- `ASMC`, written in C++, estimates pairwise coalescence times in all pairs of samples in a data set.

The `ASMCprepareDecoding` is contained in a cross-platform Jar file and does not need compiling.

To compile the `ASMC` program, edit `ASMC/SRC/Makefile` and update the location of header and binary files for the Boost library (you should not need to do this if using `getBoost.sh` was successful). Also in the Makefile, uncomment the type of [SIMD](https://en.wikipedia.org/wiki/SIMD) instruction set you want to use, depending on what is supported by your machine. Options are NO_SSE (no SIMD instructions, slower), SSE, AVX (default, recommended), AVX512. Compile the code by running "make" in the ASMC/SRC folder, or run the build.sh script from the base folder:
```
 	sh build.sh
```

### Running ASMC <a name="running" />

To run ASMC, you need to first compute a set of *decoding quantities* using the `ASMCprepareDecoding` program. Once you have obtained the decoding quantities, you can use them to analyze SNP or WGS data using the `ASMC` program.

##### Running ASMCprepareDecoding

The `ASMCprepareDecoding` precomputes all the information needed by `ASMC` for a given set of parameters (the *decoding quantities*). Detailed command line options are described [here](#detailed).

If you are analyzing European data and would like to use one of the time discretizations contained in the `FILES/DISC/*` folder, you don't need to run `ASMCprepareDecoding`. You can download the corresponding set of decoding quantities from `FILES/DECODING_QUANTITIES/*` (details [here](#precomputed)).

The following command builds decoding quantities for a European demogramic model, UK-Biobank SNP allele frequencies, and the `30-100-2000.disc` time discretization provided in the `FILES/` folder:
```
java -jar TOOLS/PREPARE_DECODING/ASMCprepareDecoding.jar \
        --demography FILES/CEU.demo \
        --discretization FILES/DISC/30-100-2000.disc \
        --freqFile FILES/UKBB.frq \
        --samples 100 \
        --out FILES/EXAMPLE/exampleFile.n100
```

The most important arguments needed by `ASMCprepareDecoding` to compute decoding quantities are as follows (file formats are described [here](#formats)):

- A file containing the *demographic history* of the analyzed samples (`--demography`). If not specified, a default European demographic model is assumed. The demographic model should roughly match that of the analyzed population, although it needs not be extemely accurate
- The set of *discrete time intervals* in which all coalescence events are assumed to occur (`--discretization`). Alternatively, you can specify the number of time intervals using the `--coalescentQuantiles` flag, and the discretization will be computed internally using quantiles of the pairwise coalescence distribution.
- *SNP allele frequencies* (`--freqFile`). These are used to deal with the non-random ascertainment of SNPs in array data. You can input precomputed allele frequencies using the `--freqFile` flag (recommended). Alternatively, you can provide the path of a raw `.haps` file using the `--fileRoot` flag. Note, however, that the goal here is to compute the genome-wide allele frequency spectrum of the SNP array data, rather than the spectrum for a specific region.
- The number of samples to be used in the CSFS. The default is n=300 and can be used in most analyses. This is specified here because the data set we will analyze contains only 100 samples.
- The root for output files (`--out`).

The program will output two files: `outFileRoot.decodingQuantities.gz` and `outFileRoot.intervalsInfo`, described [here](#formats).

##### Running ASMC

Once you have computed or downloaded the decoding quantities corresponding to your demographic model, time discretization, and SNP allele frequencies, you can analyze SNP or WGS data using the `ASMC` program:

```
ASMC/ASMC \
        --decodingQuantFile FILES/EXAMPLE/exampleFile.n100.decodingQuantities.gz \
        --hapsFileRoot FILES/EXAMPLE/exampleFile.n100.array \
        --posteriorSums
```

For WGS data, add the `--mode sequence` option:

```
ASMC/ASMC \
        --decodingQuantFile FILES/EXAMPLE/exampleFile.n100.decodingQuantities.gz \
        --hapsFileRoot FILES/EXAMPLE/exampleFile.n100 \
        --posteriorSums \
        --mode sequence
```

Using the `--majorMinorPosteriorSums` flag, ASMC will output the sum of posterior coalescence probabilities for all analyzed pairs of individuals. These will be written in `FILES/EXAMPLE/exampleFile.n100.array.{00,01,11}.sumOverPairs.gz` for the SNP array example, and `FILES/EXAMPLE/exampleFile.n100.{00,01,11}.sumOverPairs.gz` for the WGS example.

If you are decoding a large number of samples, you can break down the computation in several *jobs* using the `--jobs int` and `--jobInd int` flags, which take the total number of jobs to be performed and the current job index as arguments. If you don't speficy a name for the output files, a job index will be automatically added to the default output path. You can merge and normalize the results from all jobs using the `ASMCmergePosteriorSums` tool described [here](#mergetool).

### A complete example <a name="example" />

You can run `sh build.sh` to download and compile Boost, build ASMC (or do this [manually](#compiling)).

You can use the following commands to prepare and run ASMC using 10 jobs, and then merge the results:
```
# prepare the decoding quantities:
sh prepare.sh
# analyze array data in 10 batches:
sh decode.sh
# merge the results of the 10 batches:
sh merge.sh
```

If you want to clean things up, you could run the `sh clean.sh` script.

### Detailed command line options <a name="detailed" />

The full set of command line options for `ASMCprepareDecoding` is as follows:

```
List of arguments:
  Short            Long                         Explanation
  -h               --help                       Display help message
  -d file          --discretization file        File with time intervals
  -qCoal int       --coalescentQuantiles int    Number of generated time intervals
  -o fileRoot      --out fileRoot               Root of output files
  -D file          --demography file            File with demographic model
  -C file          --CSFS file                  File with precomputed CSFS
  -n int           --samples int                Number of samples in CSFS
  -mu float        --mut float                  Mutation rate used in demography
  -F file          --freqFile file              File with SNP allele frequencies
  -f fileRoot      --fileRoot fileRoot          Root of file with data to compute frequencies

Mandatory arguments:
  Must specify an option for time discretization:
        (-d file | -qCoal int)
  Must specify an option for SNP allele frequencies:
        (-f fileRoot | -F file)
  Must specify root of output files:
        -o fileRoot
```

In addition to the arguments described above, `ASMCprepareDecoding` options include:
- `-C` or `--CSFS`, which takes a file as argument, and allows a user to load a preloaded CSFS file (e.g. generated using TOOLS/PREPARE_DECODING/getCSFS.sh or TOOLS/PREPARE_DECODING/get_csfs.py).
- `-n` or `--samples`, which takes an integer as argument. This is used to specify how many samples to be used in the CSFS calculations. The default is 300. Any number between `100` and `300` is reasonable, with the constraint that this should be at most equal to the number of samples contained in the data set you will analyze.
- `-mu` or `--mut` is the mutation rate assumed when computing the demographic model (default=1.65E-8).

The full set of command line options for `ASMC` is as follows:

```
Mandatory:
  --hapsFileRoot arg            Prefix of hap|haps|hap.gz|haps.gz and sample|samples file
  --decodingQuantFile arg       Decoding quantities file

Choose one of:
  --posteriorSums               Output file for sum of posterior distribution
                                over pairs.
  --majorMinorPosteriorSums     Output file for sum of posterior distribution 
                                over pairs, partitioned by major/minor alleles.                                

Optional:
  --outFileRoot arg             Output file for sum of posterior distribution
                                over pairs (default: --hapsFileRoot argument)
  --jobs int (=0)               Number of jobs being done in parallel
  --jobInd int (=0)             Job index (1..jobs)
  --mode string (=array)        Decoding mode. Choose from {sequence, array}.
  --compress (=false)           Compress emission to binary (no CSFS)
  --useAncestral (=false)       Assume ancestral alleles are coded as 1 in
                                input (will assume 1 = minor otherwise)
  --skipCSFSdistance int (=0)   Genetic distance between two CSFS emissions
```

In addition to the arguments described above, `ASMC` options include:
- `--compress` is a shorthand for `--skipCSFSdistance Infinity` (see below).
- `--useAncestral` can be used to specify that a `1` in the data specifies an ancestrall allele. This will cause the CSFS to be used without folding. This is mostly not needed.
- `--skipCSFSdistance int`, which takes an integer argument specifies the minimum distance for a CSFS emission to be used. The default is `0` (always use CSFS). Setting `--skipCSFSdistance Infinity` (which is the same as `--compress`) leads to never using the CSFS (i.e. the classic PSMC emission if decoding WGS data, or a binary emission which controls for ascertainment if decoding SNP array data).

### Input/output file formats <a name="formats" />

You may want to look at files in FILES/\* for examples of the file formats described below.

##### Phased haplotypes in [Oxford haps/sample](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample) format (\*.hap/hap.gz, samples)

These files are provided in input to `ASMC` and optionally `ASMCprepareDecoding`. The file format explained [here](http://www.shapeit.fr/pages/m02_formats/hapssample.html), see "SAMPLE file" and "HAPS file" sections. These files are output by phasing programs like Eagle and Shapeit.

##### Genetic map (\*.map/map.gz)

The genetic map provided in input to `ASMC` is in [Plink map](https://www.cog-genomics.org/plink2/formats#map) format, in which each line has four columns with format "Chromosome SNPName GeneticPosition PhysicalPosition". Genetic positions are in centimorgans, physical positions are in bp. The map can be optionally compressed using gzip.

##### Demographic history (\*.demo)

The demographic history provided in input to `ASMCprepareDecoding` represents a piece-wise constant history of past effective population sizes, with format
```
TimeStart   PopulationSize
```

Where TimeStart is the first generation where the population has size PopulationSize. Note that population size is *haploid*, and that the demographic model is usually built assuming a specific mutation rate, which is passed as an argument to the `ASMCprepareDecoding` program. The first line should contain generation `0`. You can obtain this model using e.g. PSMC/MSMC/SMC++. If your model is not piecewise constant, you will need to approximate it as piecewise constant.

##### Time discretization (\*.disc)

The list of discrete time intervals provided in input to `ASMC` contains a single number per line, representing time measured in (continuous) generations, and starting at generation `0.0`. For instance, the list `FILES/DISC/10.disc` contains 10 time intervals:
```
0.0
1118.2
1472.2
1849.7
2497.0
3963.8
9120.8
15832.9
24139.9
34891.6
```

The intervals defined in this file are: `{0.0-1118.2, 1118.2-1472.2, 1472.2-1849.7, 1849.7-2497.0, 2497.0-3963.8, 3963.8-9120.8, 9120.8-15832.9, 15832.9-24139.9, 24139.9-34891.6, 34891.6-Infinity}`.

##### Decoding quantities (\*.decodingQuantities.gz)

The `*.decodingQuantities.gz` file is generated by `ASMCprepareDecoding` and input into `ASMC`. It is used to perform efficient inference of pairwise coalescence times. There is no need to understand the content of this file.

##### Time discretizzation intervals (\*.intervalsInfo)

The `*.intervalsInfo` file is generated by the `ASMCprepareDecoding` and input into `ASMC`. It contains some useful information about the time discretization and the demographic model. It contains a number of lines corresponding to the number of discrete time intervals used in the analysis. Each line has format:

```
IntervalStart   ExpectedCoalescenceTime IntervalEnd
```

IntervalStart and IntervalEnd represent the start/end of each discrete time interval, ExpectedCoalescenceTime is the expected coalescence time for a pair of individuals who have been inferred to coalesce within this time interval, and depends on the demographic model.

##### Sum of pairwise posterior coalescence probabilities `*.{00,01,11}.sumOverPairs.gz`

The output of the `ASMC` analysis is written in `*.{00,01,11}.sumOverPairs.gz` files. Each file contains a matrix of size SxD, where S is the number of sites in the data, and D is the number of discrete time intervals used in the analysis. The [*i*,*j*]-th entry of each matrix contains the sum of posterior coalescence probabilities for all samples at SNP *i* and discrete coalescence time *j*. The output breaks down coalescence events of samples carrying different alleles at each site, using the `{00,01,11}` suffixes. Specifically:

- The *i*-th row of the matrix in `*.00.sumOverPairs.gz` contains the sum of posterior probabilities for all pairs of samples that are homozygous `0` at site *i*.
- The *i*-th row of the matrix in `*.01.sumOverPairs.gz` contains the sum of posterior probabilities for all pairs of samples that heterozygous at site *i*.
- The *i*-th row of the matrix in `*.11.sumOverPairs.gz` contains the sum of posterior probabilities for all pairs of samples that are homozygous `1` at site *i*.

### Tools and scripts <a name="tools" />

These are some useful tools and scripts to be used with ASMC.

#### Tool to merge output of parallel ASMC jobs

The folder `TOOLS/MERGE_POSTERIORS/` contains the `ASMCmergePosteriorSums.jar` program, which may be used to merge the output of different ASMC jobs. You can type `java -jar TOOLS/MERGE_POSTERIORS/ASMCmergePosteriorSums.jar  -h` for a list of comman line options. Also see the `merge.sh` file. This tool assumes the decoding has been done using `--majorMinorPosteriorSums`. You may normalize the output so that the posterior sums to `1` for each site using the `--norm` flag. If you used the `--posteriorSums` flag and you want to simply normalize the output, you can run:

```
zcat FILES/EXAMPLE/exampleFile.n100.array.merged.sumOverPairs.gz | \
    awk '{ c=0; for (i=1; i<=NF; i++) c+=$i; for (i=1; i<=NF; i++) $i/=c; print; }' | \
    gzip -c -v - > FILES/EXAMPLE/exampleFile.n100.array.merged.norm.sumOverPairs.gz
```

#### Script to visualize average coalescence density in a region

This tool can be used to visualize coalescence density in specific regions (e.g. Figure 3.b/c in the ASMC paper). The script will be uploaded soon, together with results from the UKBB and GoNL analyses.

### Precomputed decoding quantities <a name="precomputed" />

[This folder](http://www.stats.ox.ac.uk/~palamara/ASMC/decodingquantitites/) contains several sets of decoding quantities that have been precomputed by running
```
sh FILES/DECODING_QUANTITIES/generate.sh
```
They are built using the European demographic model `FILES/CEU.demo`, SNP allele frequencies from the UK Biobank in `FILES/UKBB.frq`, and the time discretizations that can be found in `FILES/DISC/*.disc`. Some of these decoding quantities can be also found in the `FILES/DECODING_QUANTITIES/` folder.

The `FILES/DECODING_QUANTITIES/30-100-2000.disc` and `FILES/DECODING_QUANTITIES/10-20-2000.disc` files contain several short time intervals in the recent generations, and the same time intervals as in 60.disc from generation `2,000` on. This enables getting more fine-grained information for recent generation, though note that smaller time intervals will contain fewer coalescent events on average.

### Change log <a name="log" />

- July 1, 2018 -- Release of ASMC v1.0.

<!-- ### Regions <a name="regions" />

The following regions are delimited by chromosome start/end or centromeres, and were analyzed in the ASMC paper:

```
Chr   From            To
1     752721          121475791
1     143199864       249222527
2     11944           92191582
2     95342412        243041411
3     72365           90485813
3     93529128        197848556
4     71566           49606599
4     52684820        190906015
5     29523           46403633
5     49441779        180698413
6     203397          58746673
6     61934516        170907734
7     41892           58013387
7     61070337        159124173
8     164984          43752221
8     46936447        146292681
9     16967           47204896
9     65507526        141101939
10    93502           39074656
10    42423532        135440226
11    193146          51562420
11    54710133        134945120
12    190980          34827687
12    37858073        133831319
13    19020095        115103150
14    19264875        107287663
15    20044342        102397317
16    85629           35257261
16    46501717        90170095
17    6157            22217629
17    25278811        81103682
18    11358           15318075
18    18546896        78015180
19    110783          24545657
19    27742769        59097933
20    63231           26246166
20    29471031        62915231
21    14595742        48099610
22    16057417        51193629
``` -->

