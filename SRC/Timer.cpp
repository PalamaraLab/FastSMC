#include <cstdlib>
#include <sys/time.h>

#include "Timer.hpp"

Timer::Timer(void) {
  update_time();
}

double Timer::update_time(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  prevtime = curtime;
  curtime = tv.tv_sec + 1e-6 * tv.tv_usec;
  return curtime - prevtime;
}

unsigned long long Timer::rdtsc(void) {
  unsigned int hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ((unsigned long long) lo) | (((unsigned long long) hi) << 32);
}

double Timer::get_time(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}
