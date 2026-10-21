#include "normaldist.h"
#include "api.h"
#include "poly.h"
#include "randombytes.h"
#define _USE_MATH_DEFINES
#include <math.h>
//#include <x86intrin.h>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <pthread.h>

//pthread_t chain[NUM_CHAINS];

void write_array_to_file(int8_t f[ANTRAG_D], const char *filename) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    fprintf(file, "[");
    for (int i = 0; i < ANTRAG_D; i++) {
        fprintf(file, "%d", f[i]);
        if (i < ANTRAG_D - 1) {
            fprintf(file, ", ");
        }
    }
    fprintf(file, "]\n");

    fclose(file);
}

static void simple_frand(double *r, uint64_t *buf, size_t n) {
    // static const double pow2m64 = pow(2,-64);
    static const double pow2m64 = 5.42101086242752217e-20;
    randombytes((uint8_t*)buf, n*sizeof(uint64_t));
    for(size_t i=0; i<n; i++) {
	     r[i] = ((double)buf[i]) * pow2m64;
    }
}

// generates single random value
double single_rng() {
    uint64_t rand_int = get64();
    double rand_double = (double)rand_int / (double)UINT64_MAX;
    return rand_double;
}

/* 
 * Decode a real vector as an integer vector with odd sum.
 *
 * In other words: solve CVP in the non trivial coset of the D_n lattice.
 *
 * This is done using the decoder described in Conway & Sloane 20.2:
 * round all coefficients to the nearest integer, and if the sum is even,
 * round the worst coefficient (the one farthest from Z) in the other
 * direction.
 *
 * This implementation below is deliberately not constant time (since we
 * are not aiming for a constant-time keygen), but this is
 * straightforward to fix if deemed necessary.
 */
static void decode_odd(int8_t u[ANTRAG_D], const poly *utilde)
{
    uint8_t umod2 = 0;
    int8_t ui, wi = 0;
    size_t worst_coeff = 0;
    double maxdiff = -1, uitilde, diff;

    for(size_t i=0; i<ANTRAG_D; i++) {
    	uitilde = utilde->coeffs[i].v;
    	ui      = lrint(uitilde);
    	umod2  ^= ui;
    	diff    = fabs(uitilde - (double)ui);
    	if(diff > maxdiff) {
    	    worst_coeff = i;
    	    maxdiff = diff;
    	    if(uitilde > (double)ui)
    		    wi = ui + 1;
    	    else
    		    wi = ui - 1;
    	}
    	u[i] = ui;
    }
    if((umod2 & 1) == 0)
	u[worst_coeff] = wi;
}

static void Simple_Rounding(int8_t u[ANTRAG_D], const poly *utilde)
{
    int8_t ui;
    
    for (size_t i = 0; i < ANTRAG_D; i++) {
        ui = lrint(utilde->coeffs[i].v);
        u[i] = ui;
    }
}

double update_xi(int d, double alpha, double *x, int i, double _r) {
    double sum_other_x = 0.0;
    double sum_other_x_inv = 0.0;

    for (int j = 0; j < d / 2; j++) {
        if (j != i) {
            sum_other_x += x[j];
            sum_other_x_inv += (1.0 / x[j]);
        }
    }

    double lower_bound = 1.0 / (d / 2.0 * alpha * alpha - sum_other_x_inv);
    double upper_bound = d / 2.0 * alpha * alpha - sum_other_x;

	return lower_bound + (upper_bound - lower_bound) * _r;
}

/**
 * Function to update coefficients xi and xj based on joint 
 * conditional distribution conditioned on other coefficients
 * 
 * author: Lauren Kaestle <lk2958@rit.edu>
 */
void update_xi_xj(int d, double alpha, double *x, int i, int j) {
    double sum_other_x = 0.0;
    double sum_other_x_inv = 0.0;
    double xi, xj;

    for (int k = 0; k < d / 2; k++) {
        if (k != i && k != j) {
            sum_other_x += x[k];
            sum_other_x_inv += (1.0 / x[k]);
        }
    }

    double max_sum = (d / 2.0) * alpha * alpha;
    
    // find min and max xi and xj to build sampling triangle
    double bound1 = max_sum - sum_other_x;
    double bound2 = max_sum - sum_other_x_inv;
    double x_min = (bound1 - sqrt(bound1 * (bound1 - 4.0 / bound2))) / 2.0;
    double x_max = (bound1 + sqrt(bound1 * (bound1 - 4.0 / bound2))) / 2.0;

    while(1) {
        // sample uniformly within triangle, reject if hyperbola requirement not met
        // A = (xmin, xmin), B=(xmin, xmax), C=(xmax, xmin)
        // P = (1-sqrt(r1))A + (sqrt(r1)(1-r2))B + (r2sqrt(r1))C
        double r1 = single_rng();
        double r2 = single_rng();
        //xi = (max_sum - sum_other_x) * r2 * sqrt(r1);
        //xj = (max_sum - sum_other_x) * sqrt(r1) * (1.0 - r2);
        xi = x_min * ((1.0 - sqrt(r1)) + sqrt(r1) * (1.0 - r2)) + x_max * r2 * sqrt(r1);
        xj = x_min * ((1.0 - sqrt(r1)) + r2 * sqrt(r1)) + x_max * sqrt(r1) * (1.0 - r2);
        if (xi > 0 && xj > 0 && 1.0 / xi + 1.0 / xj <= max_sum - sum_other_x_inv) {
            x[i] = xi;
            x[j] = xj;
            break;
        }
    }
}

/**
 * Function to carry out a sequential DE-MC update
 * 
 * For each coefficient, the function computes the difference between the two 
 * randomly selected chains and uses this to compute what the update would be.
 * If the updated coefficients are all still > 0 and meet the sum requirements,
 * the update is carried out. Otherwise, there's no update for this iteration.
 * 
 * author: Lauren Kaestle <lk2958@rit.edu>
 */
void update_demc(int d, double alpha, double *x, double *a, double *b, double e) {
    double gamma = 2.38 / sqrt(d);
    double sum_x = 0.0;
    double sum_x_inv = 0.0;
    double diff = 0.0;
    double max_sum = (d / 2) * alpha * alpha;
    for (int i = 0; i < d / 2; i++) {
        diff = x[i] + gamma * (a[i] - b[i]) + e;
        if (diff <= 0.0) {
            return;
        }
        sum_x += diff;
        sum_x_inv += 1.0 / diff;
        if (sum_x > max_sum || sum_x_inv > max_sum) {
            return; // no update!
        }
    }
    for (int i = 0; i < d / 2; i++) {
        x[i] += gamma * (a[i] - b[i]) + e;
    }
}

/**
 * Function for generating the secret key using an approach other than FALCON's default
 * try-and-repeat method. Edited to check which algorithm should be used and run the
 * coefficient updating algorithm with the selected approach.
 * 
 * Edited by Lauren Kaestle <lk2958@rit.edu>
 */
int keygen_fg(secret_key *sk, enum Algorithm alg)
{
    double z[ANTRAG_D/2], af[ANTRAG_D/2], ag[ANTRAG_D/2], f[ANTRAG_D], g[ANTRAG_D];
    double ANTRAG_ALPHA_PRIME = ANTRAG_ALPHA + 0.005;
    double Exit_bound = ANTRAG_D * 0.5 * ((ANTRAG_ALPHA_PRIME)*(ANTRAG_ALPHA_PRIME)) * ((double)ANTRAG_Q);
    
    int trials=0, check;
    double x[ANTRAG_D / 2];
    for (int i = 0; i < ANTRAG_D / 2; i++) {
        x[i] = 1.0;
    }

	double _r[ANTRAG_D / 2];
	uint64_t _rint[ ANTRAG_D / 2];
	
    // edited to allow different approaches:
    if (alg == GIBBS) {
        for (int k = 0; k < Gibbs_N; k++) {
		    simple_frand(_r, _rint, ANTRAG_D / 2);
            for (int i = 0; i < ANTRAG_D / 2; i++) {
                x[i] = update_xi(ANTRAG_D, ANTRAG_ALPHA, x, i, _r[i]);
            }
        }
    } else if (alg == BLOCK) {
        seed_rng();
        for (int k = 0; k < Blocked_N; k++) {
            for (int i = 0; i < ANTRAG_D / 2; i += 2) {
                update_xi_xj(ANTRAG_D, ANTRAG_ALPHA, x, i, i+1);
            }
        }
    } else if (alg == DE) {
        seed_rng();
        // DE-MC, done sequentially with NUM_CHAINS chains (for now)
        double chains[NUM_CHAINS][ANTRAG_D / 2];
        // 1. for each chain, perform an initial update using update_xi
        for (int k = 0; k < NUM_CHAINS; k++) {
            simple_frand(_r, _rint, ANTRAG_D / 2);
            for (int i = 0; i < ANTRAG_D / 2; i++) {
                // increase or decrease each element by a little bit
                chains[k][i] = x[i] + (_r[i] * DE_START * 2) - DE_START;
            }
        }
        // 2. for DE_N iterations, use demc_update
        for (int k = 0; k < DE_N; k++) {
            for (int i = 0; i < NUM_CHAINS; i++) {
                // randomly select two other chains for the update
                int i1, i2;
                do {
                    i1 = rand() % NUM_CHAINS;
                } while (i1 == i);
                do {
                    i2 = rand() % NUM_CHAINS;
                } while (i2 == i || i2 == i1);
                // randomly generate the amount of error for this update
                double err = single_rng() * DE_ERR;
                update_demc(ANTRAG_D, ANTRAG_ALPHA, chains[i], chains[i1], chains[i2], err);
            }
        }
        // 3. select vector to use, start by taking avg for each coord
        double sum_x = 0.0;
        double sum_x_inv = 0.0;
        double max_sum = (ANTRAG_D / 2.0) * ANTRAG_ALPHA * ANTRAG_ALPHA;
        for (int i = 0; i < ANTRAG_D / 2; i++) {
            x[i] = 0.0;
            for (int j = 0; j < NUM_CHAINS; j++) {
                x[i] += chains[j][i];
            }
            x[i] /= DE_N;
            sum_x += x[i];
            sum_x_inv += 1.0 / x[i];
        }
        if (sum_x > max_sum || sum_x_inv > max_sum) {
            // for now, if this fails to meet requirements, just select chain randomly
            // going forward, could also use a metric to pick "best" chain
            int i = single_rng() * NUM_CHAINS;
            for (int j = 0; j < ANTRAG_D / 2; j++) {
                x[j] = chains[i][j];
            }
        }
    } else {
        return -1;
    }

    // For testing that coefficients generated meet requirements (will delete after testing):
    printf("TEST:\n");
    double test_sum = 0.0;
    double test_inv = 0.0;
    for (int i = 0; i < ANTRAG_D/2; i++) {
        printf("%f ", x[i]);
        test_sum += x[i];
        test_inv += 1.0 / x[i];
    }
    printf("\n%f, %f, %f\n", test_sum, test_inv, (ANTRAG_D / 2.0) * ANTRAG_ALPHA * ANTRAG_ALPHA);
    

    do {
        trials++;
        simple_frand(_r, _rint, 2*ANTRAG_D);
        // e^(iω)=cos⁡(ω)+i*sin(ω)

        for(size_t i=0; i<ANTRAG_D/2; i++) {
            z[i] = sqrt((double) ANTRAG_Q * x[i]);

            af[i]            = z [i] * cos(M_PI/2*_r[i]);              // x
            ag[i]            = z [i] * sin(M_PI/2*_r[i]);              // y
            f [i]            = af[i] * cos(2*M_PI*_r[i+1*ANTRAG_D/2]); // Real part of z
            f [i+ANTRAG_D/2] = af[i] * sin(2*M_PI*_r[i+1*ANTRAG_D/2]); // Imaginary part of z
            g [i]            = ag[i] * cos(2*M_PI*_r[i+2*ANTRAG_D/2]); // Real part of z prime
            g [i+ANTRAG_D/2] = ag[i] * sin(2*M_PI*_r[i+2*ANTRAG_D/2]); // Imaginary part of z prime
        }

        for(size_t i=0; i<ANTRAG_D; i++) {
            sk->b10.coeffs[i].v = f[i];
            sk->b11.coeffs[i].v = g[i];
        }
        invFFT(&sk->b10);
        invFFT(&sk->b11);

        // Convert real vector to integer vector
        decode_odd(sk->f, &(sk->b10));
        decode_odd(sk->g, &(sk->b11));
        //Simple_Rounding(sk->f, &(sk->b10));
        //Simple_Rounding(sk->g, &(sk->b11));

        for(size_t i=0; i<ANTRAG_D; i++) {
            sk->b10.coeffs[i].v = (double)sk->f[i];
            sk->b11.coeffs[i].v = (double)sk->g[i];
        }

        // Compute the exit condition
        FFT(&sk->b10); // f
        FFT(&sk->b11); // g
        check = 0;

        double zi_sum = 0.0;
        double zi_inv_sum = 0.0;

        for(size_t i=0; i<ANTRAG_D/2; i++) {
            double zi = 
              fpr_sqr(sk->b10.coeffs[i]).v + fpr_sqr(sk->b10.coeffs[i+ANTRAG_D/2]).v +
              fpr_sqr(sk->b11.coeffs[i]).v + fpr_sqr(sk->b11.coeffs[i+ANTRAG_D/2]).v;
            zi_sum = zi_sum + zi;

            double zi_inv = 
              (ANTRAG_Q * ANTRAG_Q) / 
              (fpr_sqr(sk->b10.coeffs[i]).v + fpr_sqr(sk->b10.coeffs[i+ANTRAG_D/2]).v +
              fpr_sqr(sk->b11.coeffs[i]).v + fpr_sqr(sk->b11.coeffs[i+ANTRAG_D/2]).v);
            zi_inv_sum = zi_inv_sum + zi_inv;
        }

        if(zi_sum > Exit_bound || zi_inv_sum > Exit_bound) {
          check = 1;
        }

    } while(check); 
	// write_array_to_file(sk->f, "output.txt");
	// write_array_to_file(sk->g, "output.txt");

    return trials;
}
