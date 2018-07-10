# uncomment only one of these lines, depending on which instructions are available on your machine among SSE, AVX, or AVX512 (use NO_SSE otherwise -- code will be slower)
# uncomment for plain code (slower)
#CFLAGS += -DNO_SSE
# uncomment for SSE
#CFLAGS += -DSSE
# uncomment for AVX
CFLAGS += -DAVX -mavx
# uncomment for AVX512
#CFLAGS += -DAVX512 -xCOMMON-AVX512

# Compiler (icc for AVX512)
# CC = icc
CC = g++

# update this to change Boost library location
BOOST_INSTALL_DIR = ../BOOST/boost_1_67_0
BOOST = $(realpath $(BOOST_INSTALL_DIR))
# update this to change Boost headers and binary location
ifneq ($(strip ${BOOST}),)
	CPATHS += -I${BOOST}/
	LPATHS += -L${BOOST}/stage/lib
	LPATHS += -Wl,-rpath,${BOOST}/stage/lib
endif

#################################################################################
#################################################################################
#################################################################################

ifeq (${debug},true)
	CFLAGS += -g
else
	CFLAGS += -O3
endif
ifeq (${prof},true)
	CFLAGS += -g -pg
	LFLAGS += -pg
endif

CFLAGS += -msse -msse2 -msse3 -std=c++11 -fopenmp -Wall
LFLAGS += -fopenmp

LLIBS = -lboost_program_options -lboost_iostreams -lz
L = ${LPATHS} ${LLIBS} -lpthread -lm

T = ../BIN/ASMC
O = Data.o DecodingParams.o DecodingQuantities.o FileUtils.o Individual.o MemoryUtils.o StringUtils.o Timer.o
OMAIN = main.o $O

.PHONY: clean

$T: ${OMAIN}
	${CC} ${LFLAGS} -o $T ${OMAIN} $L -ftree-vectorize

%.o: %.cpp
	${CC} ${CFLAGS} ${CPATHS} -o $@ -c $< -ftree-vectorize
main.o: HMM.cpp

all: $T

clean:
	rm -f *.o
	rm -f $T
