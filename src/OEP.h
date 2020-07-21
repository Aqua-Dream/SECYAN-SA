#pragma once
#include "utils.h"
#include "definitions.h"

// oblivious permutation
void oblivPermutation(int* permutedIndices, uint16_t* values, int N, uint16_t* Z1, uint16_t* Z2);

// oblivious extended permutation
void OEP(int * indices, int M, int N, uint16_t* values, uint16_t* Z1, uint16_t* Z2);