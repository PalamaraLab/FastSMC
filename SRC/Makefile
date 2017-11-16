# for AVX512
# CC = icc
# all others
CC = g++

ifeq (${debug},true)
	CFLAGS += -g
else
	CFLAGS += -O3
endif
ifeq (${prof},true)
	CFLAGS += -g -pg
	LFLAGS += -pg
endif

# uncomment only one of these lines, depending on which instructions are available on your machine among SSE, AVX, or AVX512 (use NO_SSE otherwise -- code will be slower)
#CFLAGS += -DNO_SSE
#CFLAGS += -DSSE
CFLAGS += -DAVX
#CFLAGS += -DAVX512

# for AVX512
# CFLAGS += -xCOMMON-AVX512
# CFLAGS += -std=c++11 -fopenmp -Wall

# # for AVX
CFLAGS += -mavx
CFLAGS += -std=c++11 -msse -msse2 -msse3 -fopenmp -Wall

LFLAGS += -fopenmp

BOOST_INSTALL_DIR = ${HOME}/SOFTWARE/boost_1_58_0

# add Boost include and lib paths
ifneq ($(strip ${BOOST_INSTALL_DIR}),)
	CPATHS += -I${BOOST_INSTALL_DIR}/include
	LPATHS += -L${BOOST_INSTALL_DIR}/lib
	LPATHS += -Wl,-rpath,${BOOST_INSTALL_DIR}/lib
endif

# build link line (minus flags)
LLIBS = -lboost_program_options -lboost_iostreams -lz
L = ${LPATHS} ${LLIBS} -lpthread -lm


T = ASMC
O = Data.o DecodingParams.o DecodingQuantities.o FileUtils.o Individual.o MemoryUtils.o StringUtils.o Timer.o
OMAIN = main.o $O

.PHONY: clean

$T: ${OMAIN}
	${CC} ${LFLAGS} -o $T ${OMAIN} $L

%.o: %.cpp
	${CC} ${CFLAGS} ${CPATHS} -o $@ -c $<
main.o: HMM.cpp

all: $T

clean:
	rm -f *.o
	rm -f $T
