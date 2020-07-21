#pragma once
#include <cstdint>
uint64_t poly_eval(uint64_t* coeff, uint64_t x, int size);
void interpolate(uint64_t* X, uint64_t* Y, int size, uint64_t* coeff);