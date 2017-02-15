# for AVX512
# CC = icc
# all others
# CC = g++
CC = /usr/local/Cellar/gcc49/4.9.3/bin/g++-4.9

ifeq (${debug},true)
	CFLAGS += -g
else
	CFLAGS += -O3
endif
ifeq (${prof},true)
	CFLAGS += -g -pg
	LFLAGS += -pg
endif

# for AVX512
# CFLAGS += -xCOMMON-AVX512
# CFLAGS += -std=c++11 -fopenmp -Wall

# # for AVX
CFLAGS += -mavx
CFLAGS += -std=c++11 -msse -msse2 -msse3 -fopenmp -Wall

LFLAGS += -fopenmp

BOOST_INSTALL_DIR = /Users/pier/SOFTWARE/boost_1_58_0

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
