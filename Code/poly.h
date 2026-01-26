#ifndef POLY_H
#define POLY_H

#include "param.h"
#include "fpr.h"

typedef struct{
  fpr coeffs[ANTRAG_D];
} poly;


typedef struct{
  int8_t f[ANTRAG_D];
  int8_t g[ANTRAG_D];
  int8_t F[ANTRAG_D];
  int8_t G[ANTRAG_D];
  poly b10;
  poly b11;
  poly b20;
  poly b21;
  poly GSO_b10; //~b1[0]/<~b1, ~b1>
  poly GSO_b11; //~b1[1]/<~b1, ~b1>
  poly GSO_b20; //~b2[0]/<~b2, ~b2>
  poly GSO_b21; //~b2[1]/<~b2, ~b2>
  poly beta10;
  poly beta11;
  poly beta20;
  poly beta21;
  poly sigma1;
  poly sigma2;
} secret_key;

typedef struct{
  poly h;
  uint16_t hint[ANTRAG_D];
} public_key;

typedef struct{
  poly s1;
  poly s2;
  uint8_t r[ANTRAG_K/8];
} signature;


void print_poly(const poly* p);
void FFT(poly* p);
void invFFT(poly* p);
void pointwise_mul(poly* p1, const poly* p2);
void FFT_adj(poly* p);
void FFT_mul_adj(poly* p1, const poly* p2);
void FFT_mul_selfadj(poly *p1);
void poly_div_FFT(poly* p1, const poly* p2);
  

void naive_mul(poly* c, const poly* a, const poly* b);
void set_poly(poly* p1, const poly* p2);
void poly_add(poly* p1, const poly* p2);
void poly_sub(poly* p1, const poly* p2);

void poly_recenter(poly* p);
#endif
