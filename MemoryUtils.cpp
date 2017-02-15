#include <cstdlib>
#include <iostream>

#include "MemoryUtils.hpp"
#include "Types.hpp"

void *ALIGNED_MALLOC(uint64 size) {
#ifdef USE_MKL_MALLOC
  void *p = mkl_malloc(size, MEM_ALIGNMENT);
#else
  void *p = _mm_malloc(size, MEM_ALIGNMENT);
#endif
  // TODO: change to assert() or dispense with altogether and change ALIGNED_MALLOC to macro?
  if (p == NULL) {
    std::cerr << "ERROR: Failed to allocate " << size << " bytes" << std::endl;
    exit(1);
  } else if ((uint64) p & 0xf) {
    std::cerr << "ERROR: Memory alignment of " << size << " bytes failed" << std::endl;
    exit(1);
  }
  return p;
}
