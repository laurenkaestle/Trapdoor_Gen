#include "cpucycles.h"
#include <time.h>

/**
 * asm volatile failing on my laptop, for now using 
 * __builtin_readcyclecounter but will continue to research other methods
 */
int64_t cpucycles(void)
{ // Access system counter for benchmarking
  unsigned int hi, lo;

  #if defined(__i386__) || defined(__x86_64__)
  asm volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
  return ((int64_t)lo) | (((int64_t)hi) << 32);
  #else
  return (int64_t)__builtin_readcyclecounter();
  #endif
}