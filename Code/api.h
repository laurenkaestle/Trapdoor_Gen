#ifndef API_H
#define API_H

#include <stdint.h>
#include "poly.h"

/**
 * Algorithm enum allows specification of which approach to use by passing
 * one of these values into test_gibbs and keygen_fg
 * Currently implemented: GIBBS and BLOCK
 */
enum Algorithm {
  GIBBS, BLOCK, DE
};

int keygen_fg(secret_key *sk, enum Algorithm alg);
int keygen_full(secret_key *sk, public_key *pk);
void sign(const uint8_t* m, const secret_key* sk, signature* s);
int verify(uint8_t* m, public_key* pk, signature* s);
void sampler(const secret_key* sk, const poly* c1, const poly* c2, poly* v0, poly* v1);

// to use in all keygen files
void simple_frand(double *r, uint64_t *buf, size_t n);
void decode_odd(int8_t u[ANTRAG_D], const poly *utilde);

/* Constant-time macros */
#define LSBMASK(c)      (-((c)&1))
#define CMUX(x,y,c)     (((x)&(LSBMASK(c)))^((y)&(~LSBMASK(c))))
#define CFLIP(x,c)      CMUX(x,-(x),c)
#define CZERO64(x)      ((~(x)&((x)-1))>>63)

#endif
