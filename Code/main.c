#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>  
#include <math.h>  

#include "api.h"
#include "poly.h"
#include "cpucycles.h"
#include "normaldist.h"
#include "randombytes.h"

int intcmp(const void *x, const void *y)
{
    int ix = *(int*)x, iy = *(int*)y;
    return (ix>iy) - (ix<iy);
}
int uint64cmp(const void *x, const void *y)
{
    uint64_t ix = *(uint64_t*)x, iy = *(uint64_t*)y;
    return (ix>iy) - (ix<iy);
}
int doublecmp(const void *x, const void *y)
{
    double ix = *(double*)x, iy = *(double*)y;
    return (ix>iy) - (ix<iy);
}
/*
void speed(){
  secret_key sk;
  public_key pk;
  signature s; 
 
  uint8_t m[32] = {0x46,0xb6,0xc4,0x83,0x3f,0x61,0xfa,0x3e,0xaa,0xe9,0xad,0x4a,0x68,0x8c,0xd9,0x6e,0x22,0x6d,0x93,0x3e,0xde,0xc4,0x64,0x9a,0xb2,0x18,0x45,0x2,0xad,0xf3,0xc,0x61};
  
  keygen_full(&sk, &pk);

  int iter = 10000;
  clock_t start_time = clock();
  uint64_t start = cpucycles();
  
  for(int i=0; i < iter; ++i){
    sign(m, &sk, &s);
  }
  uint64_t stop = cpucycles();
  clock_t stop_time = clock();
  
  double delta = (double)(stop - start)/iter;
  double delta_time = (double)(stop_time - start_time)/iter*1.0e6/CLOCKS_PER_SEC;
  printf("Average number of cycles per signing: %f (%f us)\n", delta, delta_time);

  start_time = clock();
  start = cpucycles();
  for(int i=0; i < iter; ++i){
    verify(m, &pk, &s);
  }
  stop = cpucycles();
  stop_time = clock();
  delta = (double)(stop - start)/iter;
  delta_time = (double)(stop_time - start_time)/iter*1.0e6/CLOCKS_PER_SEC;
  printf("Average number of cycles per verification: %f (%f us)\n", delta, delta_time);
}

void speed_keygen(){
  secret_key sk;
  public_key pk;

  int iter = 100;
  clock_t start_time = clock();
  uint64_t start = cpucycles();
  
  for(int i=0; i < iter; ++i){
    keygen_full(&sk, &pk);
  }
  uint64_t stop = cpucycles();
  clock_t stop_time = clock();
  
  double delta = (double)(stop - start)/iter;
  double delta_time = (double)(stop_time - start_time)/iter*1.0e3/CLOCKS_PER_SEC;
  printf("Average number of cycles per keygen: %f (%f ms)\n", delta, delta_time);

}
*/

/**
 * Function to test speed and cycles of trapdoor sampling with a specified algorithm
 */
int test_gibbs(enum Algorithm alg){
  srand(time(0));
  seed_rng();
  printf("Hello world, signature is %s-%u\n", (alg == GIBBS) ? "Gibbs-Falcon" : "Blocked-Gibbs-Falcon", ANTRAG_D);
  secret_key sk;
  public_key pk;
  signature s;

  // function pointer based on algorithm choice
  if (alg != GIBBS && alg != BLOCK) { // add DE once implemented
    printf("Error: Unknown algorithm type\n");
    return -1;
  }

  uint8_t m[32] = {0x46,0xb6,0xc4,0x83,0x3f,0x61,0xfa,0x3e,0xaa,0xe9,0xad,0x4a,0x68,0x8c,0xd9,0x6e,0x22,0x6d,0x93,0x3e,0xde,0xc4,0x64,0x9a,0xb2,0x18,0x45,0x2,0xad,0xf3,0xc,0x61};
  
  printf("\n* Generate initial key pair.\n");
  // keygen_full(&sk, &pk);

  int res_trials = keygen_fg(&sk, alg);
  printf("  res_trials = %d\n", res_trials);
  printf("  ...done.\n\n");
 
#define KEYGEN_TESTS 100
  printf("* Test keygen repetitions (alpha=%.3f, tests=%d, N=%d).\n\n",
    ANTRAG_ALPHA, KEYGEN_TESTS, (alg == GIBBS) ? Gibbs_N : Blocked_N);

  printf("                       min  lowq  med.  uppq   max  avg.\n");
  printf("----------------------------------------------------------\n");
  fflush(stdout);
  int trials[KEYGEN_TESTS];
  double trialsavg = 0.;
  for(int i=0; i<KEYGEN_TESTS; i++) {
    trials[i] = keygen_fg(&sk, alg);
    trialsavg += trials[i];
  }
  qsort(trials, KEYGEN_TESTS, sizeof(int), intcmp);
  trialsavg /= KEYGEN_TESTS;

  printf("keygen_fg repetitions %4d %5d %5d %5d %5d %5.2f\n",
    trials[0], trials[KEYGEN_TESTS/4], trials[KEYGEN_TESTS/2],
    trials[3*KEYGEN_TESTS/4], trials[KEYGEN_TESTS-1], trialsavg);
  fflush(stdout); 

  printf("----------------------------------------------------------\n\n");

  printf("* Benchmarking the scheme.\n\n");
  printf("                       min  lowq  med.  uppq   max  avg.\n");
  printf("----------------------------------------------------------\n");

  double keygen_cycles[KEYGEN_TESTS];
  double keygen_ms[KEYGEN_TESTS];
  double keygen_cycles_avg = 0., keygen_ms_avg = 0.;
  
  for(int i=0; i < KEYGEN_TESTS; ++i){
    clock_t start_time = clock();
    uint64_t start = cpucycles();
    // keygen_full(&sk, &pk);
    keygen_fg(&sk, alg);
    uint64_t stop = cpucycles();
    clock_t stop_time = clock();

    keygen_cycles[i]   = (double)(stop - start)/1.0e6;
    keygen_ms[i]       = (double)(stop_time - start_time)*1.0e3/CLOCKS_PER_SEC;
    keygen_cycles_avg += keygen_cycles[i];
    keygen_ms_avg     += keygen_ms[i];
  }
  keygen_cycles_avg /= KEYGEN_TESTS;
  keygen_ms_avg     /= KEYGEN_TESTS;

  qsort(keygen_cycles, KEYGEN_TESTS, sizeof(uint64_t), doublecmp);
  qsort(keygen_ms,     KEYGEN_TESTS, sizeof(double),   doublecmp);
  
  printf("keygen Mcycles       %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
    keygen_cycles[0],              keygen_cycles[KEYGEN_TESTS/4],
    keygen_cycles[KEYGEN_TESTS/2], keygen_cycles[3*KEYGEN_TESTS/4],
    keygen_cycles[KEYGEN_TESTS-1], keygen_cycles_avg);
  printf("keygen speed (ms)    %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
    keygen_ms[0],              keygen_ms[KEYGEN_TESTS/4],
    keygen_ms[KEYGEN_TESTS/2], keygen_ms[3*KEYGEN_TESTS/4],
    keygen_ms[KEYGEN_TESTS-1], keygen_ms_avg);

  printf("----------------------------------------------------------\n\n");

}

int main(){
/*


  srand(time(0));
  seed_rng();
  printf("Hello world, signature is Antrag-%u\n", ANTRAG_D);
  secret_key sk;
  public_key pk;
  signature s;

  uint8_t m[32] = {0x46,0xb6,0xc4,0x83,0x3f,0x61,0xfa,0x3e,0xaa,0xe9,0xad,0x4a,0x68,0x8c,0xd9,0x6e,0x22,0x6d,0x93,0x3e,0xde,0xc4,0x64,0x9a,0xb2,0x18,0x45,0x2,0xad,0xf3,0xc,0x61};
  
  printf("\n* Generate initial key pair.\n");
  keygen_full(&sk, &pk);
  printf("  ...done.\n\n");

// #define SIGNVERIF_TESTS 10000
//   printf("* Test correctness of the scheme.\n");
//   int correct = 0;
//   for(int i=0; i<SIGNVERIF_TESTS; i++) {
//     // printf("Coming Here2");
//     // fflush(stdout);  // 强制刷新缓冲区
//     sign(m, &sk, &s);
//     // printf("Coming Here3");
//     // fflush(stdout);  // 强制刷新缓冲区
//     correct += verify(m, &pk, &s);
//   }
//   printf("  %d/%d correct signatures. (%s).\n\n", correct, SIGNVERIF_TESTS,
// 	  (correct == SIGNVERIF_TESTS)?"ok":"ERROR!");

#define KEYGEN_TESTS 1000
  printf("* Test keygen repetitions (alpha=%.2f, xi=%.3f, tests=%d).\n\n",
	  ANTRAG_ALPHA, ANTRAG_XI, KEYGEN_TESTS);

  printf("                       min  lowq  med.  uppq   max  avg.\n");
  printf("----------------------------------------------------------\n");
  int trials[KEYGEN_TESTS];
  double trialsavg = 0.;
  for(int i=0; i<KEYGEN_TESTS; i++) {
    trials[i] = _keygen_fg(&sk);
    trialsavg+= trials[i];
  }
  qsort(trials, KEYGEN_TESTS, sizeof(int), intcmp);
  trialsavg /= KEYGEN_TESTS;
   
  printf("keygen_fg repetitions %4d %5d %5d %5d %5d %5.2f\n",
	  trials[0], trials[KEYGEN_TESTS/4], trials[KEYGEN_TESTS/2],
	  trials[3*KEYGEN_TESTS/4], trials[KEYGEN_TESTS-1], trialsavg);

  // trialsavg = 0.;
  // for(int i=0; i<KEYGEN_TESTS; i++) {
  //   trials[i] = keygen_full(&sk, &pk);
  //   trialsavg+= trials[i];
  // }
  // qsort(trials, KEYGEN_TESTS, sizeof(int), intcmp);
  // trialsavg /= KEYGEN_TESTS;
  // printf("keygen repetitions   %5d %5d %5d %5d %5d %5.2f\n",
	//   trials[0], trials[KEYGEN_TESTS/4], trials[KEYGEN_TESTS/2],
	//   trials[3*KEYGEN_TESTS/4], trials[KEYGEN_TESTS-1], trialsavg);
  printf("----------------------------------------------------------\n\n");

  printf("* Benchmarking the scheme.\n\n");
  printf("                       min  lowq  med.  uppq   max  avg.\n");
  printf("----------------------------------------------------------\n");

  double keygen_cycles[KEYGEN_TESTS];
  double keygen_ms[KEYGEN_TESTS];
  double keygen_cycles_avg = 0., keygen_ms_avg = 0.;
  
  for(int i=0; i < KEYGEN_TESTS; ++i){
    clock_t start_time = clock();
    uint64_t start = cpucycles();
    _keygen_full(&sk, &pk);
    // _keygen_fg(&sk);
    uint64_t stop = cpucycles();
    clock_t stop_time = clock();

    keygen_cycles[i]   = (double)(stop - start)/1.0e6;
    keygen_ms[i]       = (double)(stop_time - start_time)*1.0e3/CLOCKS_PER_SEC;
    keygen_cycles_avg += keygen_cycles[i];
    keygen_ms_avg     += keygen_ms[i];
  }
  keygen_cycles_avg /= KEYGEN_TESTS;
  keygen_ms_avg     /= KEYGEN_TESTS;

  qsort(keygen_cycles, KEYGEN_TESTS, sizeof(uint64_t), doublecmp);
  qsort(keygen_ms,     KEYGEN_TESTS, sizeof(double),   doublecmp);
  
  printf("keygen Mcycles       %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
	  keygen_cycles[0],              keygen_cycles[KEYGEN_TESTS/4],
	  keygen_cycles[KEYGEN_TESTS/2], keygen_cycles[3*KEYGEN_TESTS/4],
	  keygen_cycles[KEYGEN_TESTS-1], keygen_cycles_avg);
  printf("keygen speed (ms)    %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
	  keygen_ms[0],              keygen_ms[KEYGEN_TESTS/4],
	  keygen_ms[KEYGEN_TESTS/2], keygen_ms[3*KEYGEN_TESTS/4],
	  keygen_ms[KEYGEN_TESTS-1], keygen_ms_avg);


*/


  // double sign_cycles[SIGNVERIF_TESTS];
  // double sign_us[SIGNVERIF_TESTS];
  // double sign_cycles_avg = 0., sign_us_avg = 0.;
  
  // for(int i=0; i < SIGNVERIF_TESTS; ++i){
  //   clock_t start_time = clock();
  //   uint64_t start = cpucycles();
  //   sign(m, &sk, &s);
  //   uint64_t stop = cpucycles();
  //   clock_t stop_time = clock();

  //   sign_cycles[i]   = (double)(stop - start)/1.0e3;
  //   sign_us[i]       = (double)(stop_time - start_time)*1.0e6/CLOCKS_PER_SEC;
  //   sign_cycles_avg += sign_cycles[i];
  //   sign_us_avg     += sign_us[i];
  // }
  // sign_cycles_avg /= SIGNVERIF_TESTS;
  // sign_us_avg     /= SIGNVERIF_TESTS;

  // qsort(sign_cycles, SIGNVERIF_TESTS, sizeof(uint64_t), doublecmp);
  // qsort(sign_us,     SIGNVERIF_TESTS, sizeof(double),   doublecmp);
  
  // printf("sign   kcycles       %5.f %5.f %5.f %5.f %5.f %5.f\n",
	//   sign_cycles[0],                 sign_cycles[SIGNVERIF_TESTS/4],
	//   sign_cycles[SIGNVERIF_TESTS/2], sign_cycles[3*SIGNVERIF_TESTS/4],
	//   sign_cycles[SIGNVERIF_TESTS-1], sign_cycles_avg);
  // printf("sign   speed (us)    %5.f %5.f %5.f %5.f %5.f %5.f\n",
	//   sign_us[0],                 sign_us[SIGNVERIF_TESTS/4],
	//   sign_us[SIGNVERIF_TESTS/2], sign_us[3*SIGNVERIF_TESTS/4],
	//   sign_us[SIGNVERIF_TESTS-1], sign_us_avg);

  // double verif_cycles[SIGNVERIF_TESTS];
  // double verif_us[SIGNVERIF_TESTS];
  // double verif_cycles_avg = 0., verif_us_avg = 0.;
  
  // for(int i=0; i < SIGNVERIF_TESTS; ++i){
  //   clock_t start_time = clock();
  //   uint64_t start = cpucycles();
  //   verify(m, &pk, &s);
  //   uint64_t stop = cpucycles();
  //   clock_t stop_time = clock();

  //   verif_cycles[i]   = (double)(stop - start)/1.0e3;
  //   verif_us[i]       = (double)(stop_time - start_time)*1.0e6/CLOCKS_PER_SEC;
  //   verif_cycles_avg += verif_cycles[i];
  //   verif_us_avg     += verif_us[i];
  // }
  // verif_cycles_avg /= SIGNVERIF_TESTS;
  // verif_us_avg     /= SIGNVERIF_TESTS;

  // qsort(verif_cycles, SIGNVERIF_TESTS, sizeof(uint64_t), doublecmp);
  // qsort(verif_us,     SIGNVERIF_TESTS, sizeof(double),   doublecmp);
  
  // printf("verif  kcycles       %5.f %5.f %5.f %5.f %5.f %5.f\n",
	//   verif_cycles[0],                 verif_cycles[SIGNVERIF_TESTS/4],
	//   verif_cycles[SIGNVERIF_TESTS/2], verif_cycles[3*SIGNVERIF_TESTS/4],
	//   verif_cycles[SIGNVERIF_TESTS-1], verif_cycles_avg);
  // printf("verif  speed (us)    %5.f %5.f %5.f %5.f %5.f %5.f\n",
	//   verif_us[0],                 verif_us[SIGNVERIF_TESTS/4],
	//   verif_us[SIGNVERIF_TESTS/2], verif_us[3*SIGNVERIF_TESTS/4],
	//   verif_us[SIGNVERIF_TESTS-1], verif_us_avg);
  printf("----------------------------------------------------------\n\n");

  //=======================================================================

  test_gibbs(GIBBS);
  test_gibbs(BLOCK);
  
  return 0;
}
