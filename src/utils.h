#pragma once
#include <cstdint>

// Global functions that record the total communication cost.
void addComm(uint64_t comm);
uint64_t getComm();
void getComm(char* str);
void resetComm();

// Share an annotation into two parts
void shareAnnot(uint16_t x, uint16_t& x1, uint16_t& x2);

// Fast computation of floor(log2(v)) where v is a 32-bit integer
int log2_32(uint32_t value);

//void equalityTest(bool* x1, bool* x2, int length, bool& z1, bool& z2);
void equalityTest(uint64_t x1, uint64_t x2, int length, bool& z1, bool& z2);

void shareMul(uint16_t x0, uint16_t x1, uint16_t y0, uint16_t y1, uint16_t& z0, uint16_t& z1);

void oblivSwitch(uint32_t& v3, uint32_t& v4, bool bit, uint32_t& r1, uint32_t& r2);
void oblivSwitch(uint16_t& v3, uint16_t& v4, bool bit, uint16_t& r1, uint16_t& r2);

// Alice and Bob share msg0 and msg1; the selection bit b is also shared
// After OET, they learn msgb in a shared form
void oblivExtTransfer(uint16_t msg01, uint16_t msg02, uint16_t msg11, uint16_t msg12, bool b1, bool b2, uint16_t& msgb1, uint16_t& msgb2);
void oblivExtTransfer(uint32_t msg01, uint32_t msg02, uint32_t msg11, uint32_t msg12, bool b1, bool b2, uint32_t& msgb1, uint32_t& msgb2);