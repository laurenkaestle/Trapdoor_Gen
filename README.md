# A Comparison of FALCON Trapdoor Generation Techniques

This repository contains code for multiple techniques of generating the trapdoor polynomials for the FALCON digital signature scheme.

## Code

This folder contains the implementation of the different trapdoor generation approaches in C.

For this project, I have added and edited the following files:

- keygen.c has been edited to also handle blocked Gibbs
- main.c has been edited to test both implemented approaches
- other files have been edited to handle compilation issues on my laptop

## Milestone 1

The following files have been included as part of the milestone 1 submission:

- edits to keygen.c and main.c in the Code folder

## Authors

This repository is a fork from jjyydzay/Gibbs_Falcon. It contains code from FALCON's reference implementation by Thomas Pornin and implementation of the algorithm described in the paper "Generating Falcon Trapdoors via Gibbs Sampler", by Chao Sun, Thomas Espitau, Junjie Song, Jinguang Han, and Mehdi Tibouchi.
