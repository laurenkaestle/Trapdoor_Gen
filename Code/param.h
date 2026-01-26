#ifndef PARAM_H
#define PARAM_H

#define ANTRAG_D 512
#define ANTRAG_K 320
#define ANTRAG_Q 12289
#define MSG_BYTES 32
// #define R 1.32
// #define R_SQUARE 1.7424  //R^2
// #define R 1.39
// #define R_SQUARE 1.9321  //R^2

#define Gibbs_N 2500

/* sigma^2 = r^2 * alpha^2 * q */
/* gamma^2 = slack^2 * sigma^2 * 2d */
/* slack = 1.042 */
 
#if ANTRAG_D == 512
  #define ANTRAG_LOG_D 9
  #define ANTRAG_ALPHA 1.045  // 1.15
#elif ANTRAG_D == 1024
  #define ANTRAG_LOG_D 10
  #define ANTRAG_ALPHA 1.035  //1.23
#endif

// #define ANTRAG_XI (1/3.)

// #define ANTRAG_SLACK 1.042
// #define SIGMA_SQUARE (R_SQUARE * ANTRAG_ALPHA * ANTRAG_ALPHA * ANTRAG_Q)
// #define GAMMA_SQUARE (ANTRAG_SLACK * ANTRAG_SLACK * SIGMA_SQUARE * ANTRAG_D * 2)


#endif
