#include "utils.h"
#include "definitions.h"
//#include <stdlib.h>
#include "OT.h"
#include "PRNG.h"
#include <random>
uint64_t total_comm = 0;

void addComm(uint64_t comm) { total_comm += comm; }

uint64_t getComm() { return total_comm; }

void getComm(char* str)
{
	if (total_comm < 1 << 13)
		sprintf(str, "%llu bits", total_comm);
	else if (total_comm < 1 << 23)
		sprintf(str, "%.2f KB", (total_comm >> 3) / 1024.0);
	else if (total_comm < 1ull << 33)
		sprintf(str, "%.2f MB", (total_comm >> 13) / 1024.0);
	else
		sprintf(str, "%.2f GB", (total_comm >> 23) / 1024.0);
}

void resetComm() { total_comm = 0; }

void shareAnnot(uint16_t x, uint16_t& x1, uint16_t& x2)
{
	x1 = gRNG.nextUInt16();
	x2 = x - x1;
}

const int tab32[32] = {
	 0,  9,  1, 10, 13, 21,  2, 29,
	11, 14, 16, 18, 22, 25,  3, 30,
	 8, 12, 20, 28, 15, 17, 24,  7,
	19, 27, 23,  6, 26,  5,  4, 31 };

// fast and cross-platform implementation
int log2_32(uint32_t value)
{
	value |= value >> 1;
	value |= value >> 2;
	value |= value >> 4;
	value |= value >> 8;
	value |= value >> 16;
	return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

// More Efficient Oblivious Transfer
inline void GMW_aubv(bool& a, bool& u, bool& b, bool& v)
{
	a = gRNG.nextBit();
	OT::GMW_OT(a, v, b, u);
	b ^= v;
}

void genGMW_Triple(bool& a0, bool& a1, bool& b0, bool& b1, bool& c0, bool& c1)
{
	bool u0 = 0, v0 = 0, u1 = 0, v1 = 0;
	GMW_aubv(a0, u0, b1, v1);
	GMW_aubv(a1, u1, b0, v0);
	c0 = (a0 && b0) ^ u0 ^ v0;
	c1 = (a1 && b1) ^ u1 ^ v1;
}

void GMW_Mul(bool x0, bool x1, bool y0, bool y1, bool& z0, bool& z1)
{
	bool a0 = 0, a1 = 0, b0 = 0, b1 = 0, c0 = 0, c1 = 0;
	genGMW_Triple(a0, a1, b0, b1, c0, c1);
	bool e0 = a0 ^ x0;
	bool e1 = a1 ^ x1;
	bool f0 = b0 ^ y0;
	bool f1 = b1 ^ y1;
	bool e = e0 ^ e1;
	bool f = f0 ^ f1;
	addComm(4); // Reveal e and f
	z0 = (f && a0) ^ (e && b0) ^ c0;
	z1 = (e && f) ^ (f && a1) ^ (e && b1) ^ c1;
}

// Attention: It will change the values of x1 and x2
void equalityTest(bool* x1, bool* x2, int length, bool& z1, bool& z2)
{
	if (length >= 256 || length <= 0)
		throw "Unsupported length for equality test!";
	while (length > 4) //New Protocols for Secure Equality Test and Comparison
	{
		uint32_t agg1 = 0, agg2 = 0;
		uint8_t modulus = (uint8_t)(length + 1);
		// length < 256, loglength < 9, msglen <= 1
		uint8_t msg0;
		uint8_t msg1;
		uint8_t msgb;
		std::uniform_int_distribution<> dist(0, modulus - 1);
		for (int i = 0; i < length; i++)
		{
			uint8_t a = dist(gRNG.rng);
			agg1 += a;
			msg0 = a + x1[i];
			if (msg0 >= modulus)
				msg0 -= modulus;
			msg1 = a + 1 - x1[i];
			if (msg1 >= modulus)
				msg1 -= modulus;
			OT::OT(&msg0, &msg1, 1, x2[i], &msgb);
			agg2 += msgb;
		}
		agg1 %= modulus;
		agg2 %= modulus;
		length = log2_32(length) + 1; // equal to ceil(log2(length+1))
		for (int i = 0; i < length; i++)
		{
			x1[i] = agg1 & 1;
			agg1 >>= 1;
			x2[i] = agg2 & 1;
			agg2 >>= 1;
		}
	}
	for (int i = 0; i < length; i++)
		x1[i] = !x1[i];
	while (length > 1)
	{
		int halfLength = length / 2;
		for (int i = 0; i < halfLength; i++)
			GMW_Mul(x1[2 * i], x2[2 * i], x1[2 * i + 1], x2[2 * i + 1], x1[i], x2[i]);
		if (length & 1)
		{
			x1[halfLength] = x1[length - 1];
			x2[halfLength] = x2[length - 1];
			length = halfLength + 1;
		}
		else
			length = halfLength;
	}
	z1 = x1[0];
	z2 = x2[0];

}

// Compare lower length-bits of x1 and x2
void equalityTest(uint64_t x1, uint64_t x2, int length, bool& z1, bool& z2)
{
	bool x1_arr[64], x2_arr[64];
	for (int i = 0; i < length; i++)
	{
		x1_arr[i] = x1 & 1;
		x1 >>= 1;
		x2_arr[i] = x2 & 1;
		x2 >>= 1;
	}
	equalityTest(x1_arr, x2_arr, length, z1, z2);
}

void preMulShare(uint16_t a, uint16_t b, uint16_t* u0, uint16_t* u1) {
	for (int i = 0; i < ell; ++i)
	{
		OT::shareMulCOT((b & 1), a, u0[i], u1[i]);
		b >>= 1;
		a <<= 1;
	}
}

void preMultiTuple(uint16_t& a0, uint16_t& a1, uint16_t& b0, uint16_t& b1, uint16_t& c0, uint16_t& c1) {
	a0 = gRNG.nextUInt16(); b0 = gRNG.nextUInt16();
	a1 = gRNG.nextUInt16(); b1 = gRNG.nextUInt16();
	uint16_t u0[ell], u1[ell], v0[ell], v1[ell];

	preMulShare(a0, b1, u0, u1);
	preMulShare(a1, b0, v0, v1);

	c0 = a0 * b0; c1 = a1 * b1;
	for (int i = 0; i < ell; ++i)
	{
		c0 += u0[i] + v0[i];
		c1 += u1[i] + v1[i];
	}
}

void shareMul(uint16_t x0, uint16_t x1, uint16_t y0, uint16_t y1, uint16_t& z0, uint16_t& z1)
{
	uint16_t a0, a1, b0, b1, c0, c1;
	preMultiTuple(a0, a1, b0, b1, c0, c1);

	uint16_t e0, e1, f0, f1, e, f;
	e0 = x0 - a0; f0 = y0 - b0;
	e1 = x1 - a1; f1 = y1 - b1;
	// Pi transmit its ei, fi to P(1-i)
	addComm(16 * 4);

	e = e0 + e1; f = f0 + f1; // each party both hold e & f
	z0 = f * a0 + e * b0 + c0;
	z1 = e * f + f * a1 + e * b1 + c1;
}

void oblivSwitch(uint32_t& v3, uint32_t& v4, bool bit, uint32_t& r1, uint32_t& r2) {
	// Party 1 hold v3, v4, bit
	// Party 2 hold r1, r2
	// If bit == true, two party swap their value

	// Party 2 generates two random value, and send to Party 1
	uint32_t r3, r4;
	r3 = gRNG.nextUInt32(); r4 = gRNG.nextUInt32();

	uint32_t msg0[2] = { r1 - r3, r2 - r4 };
	uint32_t msg1[2] = { r2 - r3, r1 - r4 };
	uint32_t msgb[2];
	OT::OT((uint8_t*)msg0, (uint8_t*)msg1, 8, bit, (uint8_t*)msgb);

	// Party 1 modify its value
	if (bit)
	{
		uint32_t t1, t2;
		t1 = v4 + msgb[0];
		t2 = v3 + msgb[1];
		v3 = t1; v4 = t2;

	}
	else
	{
		uint32_t t1, t2;
		t1 = v3 + msgb[0];
		t2 = v4 + msgb[1];
		v3 = t1; v4 = t2;
	}
	// Party 2 modify its value
	r1 = r3; r2 = r4;
}

void oblivSwitch(uint16_t& v3, uint16_t& v4, bool bit, uint16_t& r1, uint16_t& r2) {
	// Party 1 hold v3, v4, bit
	// Party 2 hold r1, r2
	// If bit == true, two party swap their value

	// Party 2 generates two random value, and send to Party 1
	uint16_t r3, r4;
	r3 = gRNG.nextUInt16(); r4 = gRNG.nextUInt16();

	uint16_t msg0[2] = { (uint16_t)(r1 - r3), (uint16_t)(r2 - r4) };
	uint16_t msg1[2] = { (uint16_t)(r2 - r3), (uint16_t)(r1 - r4) };
	uint16_t msgb[2];
	OT::OT((uint8_t*)msg0, (uint8_t*)msg1, 4, bit, (uint8_t*)msgb);

	// Party 1 modify its value
	if (bit)
	{
		uint16_t t1, t2;
		t1 = v4 + msgb[0];
		t2 = v3 + msgb[1];
		v3 = t1; v4 = t2;
	}
	else
	{
		uint16_t t1, t2;
		t1 = v3 + msgb[0];
		t2 = v4 + msgb[1];
		v3 = t1; v4 = t2;
	}
	// Party 2 modify its value
	r1 = r3; r2 = r4;
}

void oblivExtTransfer(uint16_t msg01, uint16_t msg02, uint16_t msg11, uint16_t msg12, bool b1, bool b2, uint16_t& msgb1, uint16_t& msgb2)
{
	// msg_(i,j) means the i-th share message hold by Pj
	oblivSwitch(msg01, msg11, b1, msg02, msg12);
	oblivSwitch(msg02, msg12, b2, msg01, msg11);
	msgb1 = msg01;
	msgb2 = msg02;
}

void oblivExtTransfer(uint32_t msg01, uint32_t msg02, uint32_t msg11, uint32_t msg12, bool b1, bool b2, uint32_t& msgb1, uint32_t& msgb2)
{
	// msg_(i,j) means the i-th share message hold by Pj
	oblivSwitch(msg01, msg11, b1, msg02, msg12);
	oblivSwitch(msg02, msg12, b2, msg01, msg11);
	msgb1 = msg01;
	msgb2 = msg02;
}