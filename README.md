## ASMC
Ascertained Sequentially Markovian Coalescent. Fast estimation of coalescence time along the genome for pairs of samples, using array or sequencing data.

### Dependencies
- install [SMC++](https://github.com/popgenmethods/smcpp) (used to generate CSFS).

### Compiling
After adjusting any parameters in the SRC/Makefile (e.g. Boost library path, choice of available SIMD instructions like SSE or AVX), run:

 	sh build.sh

I successfully compiled with AVX instruction set and GCC-4.8.5.

### Preparing the data
ASMC requires three data files in input: the phased haplotypes in Oxford HAPS/SAMPLE format (output by Eagle and Shapeit, explained [here](http://www.shapeit.fr/pages/m02_formats/hapssample.html) in the "SAMPLE file" and "HAPS file" sections), and a genetic map in Plink "map" format. Specifically, assuming the file root for all files is "file", ASMC will need
- file.hap.gz (the phased haplotypes in Oxford HAPS/SAMPLE format, optionally compressed using gzip)
- file.samples (the list of samples in Oxford HAPS/SAMPLE format)
- file.map.gz (the map file, in Plink "map" format, where each line is a variant with format "Chromosome Name GeneticPositionInCM PhysicalPositionInBp", optionally compressed using gzip)

The EXAMPLE folder contains sample data for a realistic simulation of 30 Mb from the beginning of chromosome 2 for 100 haploid samples from a European demographic model.

### Running the analysis
When running ASMC, three steps will be performed. These are all coded up in the *prepare.sh* and *decode.sh* scripts. The steps are:
1) Compute the CSFS using SMC++.
2) Compute a set of *decoding quantities* using the CSFS, a demographic model for the population, a list of time discretization intervals (see below for formats), and a file containing the frequencies of the variants (if analyzing array data).
3) Run the main analysis on the data.

The demographic model should roughly match that of the analyzed population, although it needs not be extemely accurate. If you have a demographic model for the population obtained using PSMC/MSMC/SMC++, you can write it using the format in FILES/CEU.demo. This is piece-wise constant (each line has format "StartGeneration SizeUntilNextChange"). If your model is not piecewise constant, you will need to approximate it as piecewise constant.

Lists of time discretization intervals (one generation value per line) are precomputed in the folder FILES/DISC. The *prepare.sh* script uses the discretization 30-100-2000.

Allele frequencies (used if you analyze array data) are currently contained in the file FILES/chrAll.frq (from the UK Biobank). This is in Plink --freq format. It should be computed using the whole array data set. It is used to fix SNP ascertainment biases.

To run the first two preparation steps, you can type:
```
sh prepare.sh EXAMPLE/exampleFile.array
```
The data analysis (step 3) can be divided into batches, which are run in parellel. I recommend dividing the genome in regions that exclude centromeres (see list at the bottom). Each region is to be treated independently for the analysis. To analyze the first of 10 batches in the example region, you can type:
```
sh decode.sh EXAMPLE/exampleFile.array array 1 10
```
Once you have run ASMC on all batches (e.g. 10 out of 10), you can merge the final results using
```
sh merge.sh EXAMPLE/exampleFile.array 10
```
If you want to analyze the example's sequencing data, you can run
```
sh prepare.sh EXAMPLE/exampleFile
sh decode.sh EXAMPLE/exampleFile sequence 1 10
```
(you could actually skip the first command and rename the files containing the precomputed quantities, since they will work for either sequencing or array analysis).

### Output format
The output will contain one row per variant, one colum per discrete time interval (in the example, these are defined in FILES/DISC/30-100-2000.disc). Each number is the sum of the inferred coalescence rate for a pair of chromosomes, for all analyzed pairs of chromosomes. Thus, each line will sum (roughly), to the total number of analyzed pairs (since the coalescence distribution for each pair will sum to 1). You can re-normalize each row if you like. For instance, you could do:
```
zcat EXAMPLE/exampleFile.array.merged.sumOverPairs.gz | awk '{c=0; for (i=1; i<=NF; i++) c+=$i; for (i=1; i<=NF; i++) $i/=c; print; }' | gzip -c -v - > EXAMPLE/exampleFile.array.merged.norm.sumOverPairs.gz
```

### Cleaning
There's a script to clean up, which you can run using:
```
sh clean.sh
```

### ASMC paper
P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. Under review.

### Regions
| Chr | From | To |
|:---:|:------:|:-----:|
| 1 | 752721 | 121475791 |
| 1 | 143199864 | 249222527 |
| 2 | 11944 | 92191582 |
| 2 | 95342412 | 243041411 |
| 3 | 72365 | 90485813 |
| 3 | 93529128 | 197848556 |
| 4 | 71566 | 49606599 |
| 4 | 52684820 | 190906015 |
| 5 | 29523 | 46403633 |
| 5 | 49441779 | 180698413 |
| 6 | 203397 | 58746673 |
| 6 | 61934516 | 170907734 |
| 7 | 41892 | 58013387 |
| 7 | 61070337 | 159124173 |
| 8 | 164984 | 43752221 |
| 8 | 46936447 | 146292681 |
| 9 | 16967 | 47204896 |
| 9 | 65507526 | 141101939 |
| 10 | 93502 | 39074656 |
| 10 | 42423532 | 135440226 |
| 11 | 193146 | 51562420 |
| 11 | 54710133 | 134945120 |
| 12 | 190980 | 34827687 |
| 12 | 37858073 | 133831319 |
| 13 | 19020095 | 115103150 |
| 14 | 19264875 | 107287663 |
| 15 | 20044342 | 102397317 |
| 16 | 85629 | 35257261 |
| 16 | 46501717 | 90170095 |
| 17 | 6157 | 22217629 |
| 17 | 25278811 | 81103682 |
| 18 | 11358 | 15318075 |
| 18 | 18546896 | 78015180 |
| 19 | 110783 | 24545657 |
| 19 | 27742769 | 59097933 |
| 20 | 63231 | 26246166 |
| 20 | 29471031 | 62915231 |
| 21 | 14595742 | 48099610 |
| 22 | 16057417 | 51193629 |
