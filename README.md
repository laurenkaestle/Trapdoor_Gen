# A Comparison of FALCON Trapdoor Generation Techniques

This repository contains code for multiple techniques of generating the trapdoor polynomials for the FALCON digital signature scheme.

## Code

This folder contains the implementation of the different trapdoor generation approaches in C.

For this project, I have added and edited the following files:

- keygen.c has been edited to also handle blocked Gibbs and DE-MC
- main.c has been edited to test all samplers
- param.h edited to add necessary macros for new algorithms
- other files have been edited to handle compilation issues on my laptop
- ConvergenceR.ipynb added for diagnostics to check when samplers are expected to reach convergence

## Convergence Calculations

The code for the convergence calculations that have been done so far is written in a combination of Python and R so that R's BOA package can be used. The Colab document with this code has been included in this repository.

## Authors

This repository is a fork from jjyydzay/Gibbs_Falcon. It contains code from FALCON's reference implementation by Thomas Pornin and implementation of the algorithm described in the paper "Generating Falcon Trapdoors via Gibbs Sampler", by Chao Sun, Thomas Espitau, Junjie Song, Jinguang Han, and Mehdi Tibouchi.
