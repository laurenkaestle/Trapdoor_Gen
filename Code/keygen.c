#include "normaldist.h"
#include "api.h"
#include "poly.h"
#include "randombytes.h"
#include <math.h>
#include <x86intrin.h>
#include <string.h>
#include <stdio.h>


void write_array_to_file(int8_t f[ANTRAG_D], const char *filename) {
    FILE *file = fopen(filename, "a");  // 使用"a"模式以追加的方式打开文件
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    // 写入Python列表格式
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
    static const double pow2m64 = pow(2,-64);
    randombytes((uint8_t*)buf, n*sizeof(uint64_t));
    for(size_t i=0; i<n; i++) {
	     r[i] = ((double)buf[i]) * pow2m64;
    }
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

// A random number generation function that returns a random floating point number in the range [min, max]

// double random_double(double min, double max) {
	// // printf("random_double: %lf\n", min + (max - min) * ((double) rand() / (double) RAND_MAX));
    // return min + (max - min) * ((double) rand() / (double) RAND_MAX);
// }

// 更新 X 向量中的第 i 个元素，依据给定条件分布
double update_xi(int d, double alpha, double *x, int i, double _r) {  //double?
    double sum_other_x = 0.0;
    double sum_other_x_inv = 0.0;

    // 计算 sum_other_x = sum(x[j]) for j != i
    for (int j = 0; j < d / 2; j++) {
        if (j != i) {
            sum_other_x += x[j];
            sum_other_x_inv += (1.0 / x[j]);
        }
    }

    // 计算条件分布的上下界
    double lower_bound = 1.0 / (d / 2.0 * alpha * alpha - sum_other_x_inv);
    double upper_bound = d / 2.0 * alpha * alpha - sum_other_x;
	
    // 在这个区间内均匀采样一个值
    // printf("_r: %lf\n", _r);
	// return random_double(lower_bound, upper_bound);
	return lower_bound + (upper_bound - lower_bound) * _r;
}

int keygen_fg(secret_key *sk)
{
    // double theta[ANTRAG_D/2], gamma_f[ANTRAG_D/2], gamma_g[ANTRAG_D/2]
    double z[ANTRAG_D/2], af[ANTRAG_D/2], ag[ANTRAG_D/2], f[ANTRAG_D], g[ANTRAG_D];
    double ANTRAG_ALPHA_PRIME = ANTRAG_ALPHA + 0.005;
    double Exit_bound = ANTRAG_D * 0.5 * ((ANTRAG_ALPHA_PRIME)*(ANTRAG_ALPHA_PRIME)) * ((double)ANTRAG_Q);
    // printf("ANTRAG_D is: %d, ANTRAG_Q is: %d, ANTRAG_ALPHA_PRIME is: %.3f, ANTRAG_ALPHA is: %.3f\n", ANTRAG_D, ANTRAG_Q, ANTRAG_ALPHA_PRIME, ANTRAG_ALPHA);

    double r[2*ANTRAG_D];
    uint64_t rint[2*ANTRAG_D];

    int trials=0, check;
    //---------------------- 
    // 初始化 x 向量 (d/2 个元素，初始为 1)
    double x[ANTRAG_D / 2];
    for (int i = 0; i < ANTRAG_D / 2; i++) {
        x[i] = 1.0;
    }

	// double _r[Gibbs_N * ANTRAG_D / 2];

    // uint64_t _rint[Gibbs_N * ANTRAG_D / 2];
	// simple_frand(_r, _rint, Gibbs_N * ANTRAG_D / 2);
	
    // 执行算法
	double _r[ANTRAG_D / 2];
	uint64_t _rint[ ANTRAG_D / 2];
	
    for (int k = 0; k < Gibbs_N; k++) {
		simple_frand(_r, _rint, ANTRAG_D / 2);
		// printf("Every _r: %lf\n", _r[0]);
        for (int i = 0; i < ANTRAG_D / 2; i++) {
            // 更新 x 中的第 i 个元素
            x[i] = update_xi(ANTRAG_D, ANTRAG_ALPHA, x, i, _r[i]);
        }
    }

    // 打印最终的 x 向量
    // printf("x: ");
    // for (int i = 0; i < d / 2; i++) {
    //     printf("%f ", x[i]);
    // }
    // printf("\n"); 

    //---------------------- 
    do {
        trials++;
        simple_frand(r, rint, 2*ANTRAG_D);
        // e^(iω)=cos⁡(ω)+i*sin(ω)

        for(size_t i=0; i<ANTRAG_D/2; i++) {
            z[i] = sqrt((double) ANTRAG_Q * x[i]);
            // printf("z[i]: %lf\n", z[i]);
            // printf("x[i]: %lf\n", x[i]);

            af[i]            = z [i] * cos(M_PI/2*r[i]);              // x
            ag[i]            = z [i] * sin(M_PI/2*r[i]);              // y
            f [i]            = af[i] * cos(2*M_PI*r[i+1*ANTRAG_D/2]); // Real part of z
            f [i+ANTRAG_D/2] = af[i] * sin(2*M_PI*r[i+1*ANTRAG_D/2]); // Imaginary part of z
            g [i]            = ag[i] * cos(2*M_PI*r[i+2*ANTRAG_D/2]); // Real part of z prime
            g [i+ANTRAG_D/2] = ag[i] * sin(2*M_PI*r[i+2*ANTRAG_D/2]); // Imaginary part of z prime
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
        // Simple_Rounding(sk->f, &(sk->b10));
        // Simple_Rounding(sk->g, &(sk->b11));

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
		// printf("zi_sum is %lf ", zi_sum);
		// printf("zi_inv_sum is %lf ", zi_inv_sum);
	

        if(zi_sum > Exit_bound || zi_inv_sum > Exit_bound) {
          check = 1;
        }
        
    } while(check); 
	// write_array_to_file(sk->f, "output.txt");
	// write_array_to_file(sk->g, "output.txt");

    return trials;
}
/*
int _keygen_fg(secret_key *sk)
{
    double z[ANTRAG_D/2], af[ANTRAG_D/2], ag[ANTRAG_D/2],
        f[ANTRAG_D], g[ANTRAG_D];
#if 0
    const double alow  = 1/ANTRAG_ALPHA + ANTRAG_ALPHAEPS,
           ahigh = ANTRAG_ALPHA - ANTRAG_ALPHAEPS,
#else
    const double alow  = 0.5*(ANTRAG_ALPHA + 1/ANTRAG_ALPHA) - 0.5*ANTRAG_XI*(ANTRAG_ALPHA-1/ANTRAG_ALPHA),
           ahigh = 0.5*(ANTRAG_ALPHA + 1/ANTRAG_ALPHA) + 0.5*ANTRAG_XI*(ANTRAG_ALPHA-1/ANTRAG_ALPHA),
#endif
     qlow  = ((double)ANTRAG_Q)*alow*alow,
           qhigh = ((double)ANTRAG_Q)*ahigh*ahigh,
           qlow2 = ((double)ANTRAG_Q)/(ANTRAG_ALPHA*ANTRAG_ALPHA),
           qhigh2= ((double)ANTRAG_Q)*ANTRAG_ALPHA*ANTRAG_ALPHA;

    double r[2*ANTRAG_D];
    uint64_t rint[2*ANTRAG_D];

    int trials=0, check;

    do {
        trials++;
        simple_frand(r, rint, 2*ANTRAG_D);
        // e^(iω)=cos⁡(ω)+i*sin(ω)
        for(size_t i=0; i<ANTRAG_D/2; i++) {
            z[i] = sqrt(qlow + (qhigh - qlow)*r[i]);

            af[i]            = z [i] * cos(M_PI/2*r[i+  ANTRAG_D/2]);
            ag[i]            = z [i] * sin(M_PI/2*r[i+  ANTRAG_D/2]);
            f [i]            = af[i] * cos(2*M_PI*r[i+2*ANTRAG_D/2]);
            f [i+ANTRAG_D/2] = af[i] * sin(2*M_PI*r[i+2*ANTRAG_D/2]);
            g [i]            = ag[i] * cos(2*M_PI*r[i+3*ANTRAG_D/2]);
            g [i+ANTRAG_D/2] = ag[i] * sin(2*M_PI*r[i+3*ANTRAG_D/2]);
        }
        for(size_t i=0; i<ANTRAG_D; i++) {
            sk->b10.coeffs[i].v = f[i];
            sk->b11.coeffs[i].v = g[i];
        }
        invFFT(&sk->b10);
        invFFT(&sk->b11);

        decode_odd(sk->f, &(sk->b10));
        decode_odd(sk->g, &(sk->b11));

        for(size_t i=0; i<ANTRAG_D; i++) {
            sk->b10.coeffs[i].v = (double)sk->f[i];
            sk->b11.coeffs[i].v = (double)sk->g[i];
        }

        FFT(&sk->b10);
        FFT(&sk->b11);
        check = 0;
        for(size_t i=0; i<ANTRAG_D/2; i++) {
            double zi =
                fpr_sqr(sk->b10.coeffs[i]).v + fpr_sqr(sk->b10.coeffs[i+ANTRAG_D/2]).v +
                fpr_sqr(sk->b11.coeffs[i]).v + fpr_sqr(sk->b11.coeffs[i+ANTRAG_D/2]).v;
            if(zi < qlow2 || zi > qhigh2) {
                check = 1;
                break;
            }
        }
    } while(check); 
    return trials;
}

int keygen_full(secret_key *sk, public_key *pk)
{
    int trials = 0, lim;
    uint8_t tmp[32*ANTRAG_D];
    
    lim = (1 << (Zf(max_FG_bits)[ANTRAG_LOG_D] - 1)) - 1;
	// printf("lim: %d", lim);
    while(1) {
      // printf("trials: %d\n", trials);fflush(stdout);
    	trials += keygen_fg(sk);
    	if (!Zf(compute_public)(pk->hint, sk->f, sk->g, ANTRAG_LOG_D, tmp)){
        // printf("Failed in compute_public\n");
        continue;
      }

    	if (!Zf(solve_NTRU)(ANTRAG_LOG_D, sk->F, sk->G, sk->f, sk->g, lim, (uint32_t*)tmp)){
        // printf("Failed in solve_NTRU\n");
        continue;
      }

    	break;
    }
    // printf("Finished generating f, g, F, G...\n");fflush(stdout);

    for(int i=0; i<ANTRAG_D; i++) {
    	sk->b20.coeffs[i].v = (double)sk->F[i];
    	sk->b21.coeffs[i].v = (double)sk->G[i];
    }
	// printf("F and G:\n");
	// for(int i=0; i<ANTRAG_D; i++) {
    	// printf("F[%d] is %lf  ", i, (double)sk->F[i]);
    	// printf("G[%d] is %lf\n", i, (double)sk->G[i]);
    // }
	// printf("\n");
    FFT(&sk->b20);
    FFT(&sk->b21);

    for(int i=0; i<ANTRAG_D; i++) {
      pk->h.coeffs[i].v = (double)pk->hint[i]; 
    }
    FFT(&pk->h);

    // compute_GSO(sk);
    // compute_sigma(sk);
    // compute_beta_hat(sk);

    // print_secret_key(sk);
	
	// FILE *file = fopen("output.txt", "a");
	// if (file == NULL) {
        // perror("Error opening file");
        // return;
    // }
	// fprintf(file, "=======================================");
    // fclose(file);
    return trials;
}


int _keygen_full(secret_key *sk, public_key *pk)
{
    int trials = 0, lim;
    uint8_t tmp[32*ANTRAG_D];

    lim = (1 << (Zf(max_FG_bits)[ANTRAG_LOG_D] - 1)) - 1;
    while(1) {
      trials += _keygen_fg(sk);
      if (!Zf(compute_public)(pk->hint, sk->f, sk->g, ANTRAG_LOG_D, tmp))
          continue;

      if (!Zf(solve_NTRU)(ANTRAG_LOG_D, sk->F, sk->G, sk->f, sk->g, lim, (uint32_t*)tmp))
          continue;

      break;
    }

    for(int i=0; i<ANTRAG_D; i++) {
      sk->b20.coeffs[i].v = (double)sk->F[i];
      sk->b21.coeffs[i].v = (double)sk->G[i];
    }
    FFT(&sk->b20);
    FFT(&sk->b21);

    for(int i=0; i<ANTRAG_D; i++) {
      pk->h.coeffs[i].v = (double)pk->hint[i]; 
    }
    FFT(&pk->h);

    // compute_GSO(sk);
    // compute_sigma(sk);
    // compute_beta_hat(sk);

    // print_secret_key(sk);
    return trials;
}
*/
// N = 200, 1000 times, chi-square test
// uniform 