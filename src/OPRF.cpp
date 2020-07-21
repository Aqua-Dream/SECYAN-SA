#include "OPRF.h"
#include "PRNG.h"
#include "aes.h"
#include "utils.h"
#include <bitset>
#include "OT.h"
#include "blake3.h"
#include <iostream>

using namespace std;
const int OPRF_k = 448;

void bitsetToArr(bitset<OPRF_k>& bs, uint8_t *arr)
{
	for (int i = 0; i < OPRF_k / 8; i++)
	{
		uint8_t v = 0;
		for (int j = 0; j < 8; j++)
		{
			v <<= 1;
			v |= bs[8 * i + j];
		}
		arr[i] = v;
	}
}

void arrToBitset(uint8_t *&arr, bitset<OPRF_k>& bs)
{
	for (int i = 0; i < OPRF_k / 8; i++)
	{
		bs <<= 8;
		bs |= arr[i];
	}
}

bitset<OPRF_k> arrToBitset(uint8_t *arr)
{
	bitset<OPRF_k> bs;
	arrToBitset(arr, bs);
	return bs;
}

bitset<OPRF_k> PRC(int index, uint64_t value)
{
	uint8_t aes_in[16];
	uint8_t aes_out[16];
	((uint64_t*)aes_in)[1] = value;
	((uint32_t*)aes_in)[1] = index;
	bitset<OPRF_k> out;
	for (int i = 0; i < 4; i++)
	{
		((uint32_t*)aes_in)[0] = i;
		aes128_enc(aes_in, aes_out);
		for (int j = 0; j < 448 / 4 / 8; j++)
		{
			out <<= 8;
			out |= aes_out[j];
		}
	}
	return out;
}

uint64_t OPRF_H(bitset<OPRF_k>& bs)
{
	const int len = OPRF_k / 8;
	uint8_t in[len];
	uint64_t out;
	bitsetToArr(bs, in);
	blake3_hasher hasher;
	blake3_hasher_init(&hasher);
	blake3_hasher_update(&hasher, in, len);
	blake3_hasher_finalize(&hasher, (uint8_t*)&out, 64/8);
	return out;
}


// We currently assume the inputs and outputs are all 64-bit length

void batch_OPRF(vector<uint32_t> &eleX, vector<vector<uint32_t>>& elesY, vector<uint64_t> &outX, vector<vector<uint64_t>>& outsY)
{
	// The sender holds elesY while the receiver holds eleX
	int m = (int)eleX.size();
	uint8_t aes_key[16];
	
	bitset<OPRF_k> s;
	for (int i = 0; i < 4; i++)
		((uint32_t*)aes_key)[i] = gRNG.nextUInt32();
	aes128_load_key_enc_only((uint8_t*)aes_key);
	addComm(128); // Sends aes_key

	// Sender chooses s at random
	for (int i = 0; i < OPRF_k / 32; i++)
	{
		s <<= 32;
		s |= gRNG.nextUInt32();
	} 
	bitset<OPRF_k>* t0 = new bitset<OPRF_k>[m];
	bitset<OPRF_k>* t1 = new bitset<OPRF_k>[m];
	bitset<OPRF_k>* q = new bitset<OPRF_k>[m];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < OPRF_k / 32; j++)
		{
			t0[i] <<= 32;
			t0[i] |= gRNG.nextUInt32();
		}
		t1[i] = t0[i] ^ PRC(i, eleX[i]);
	}
	
	int msglen = (m + 7) / 8;
	int extbits = m & 7;
	uint8_t* msg0 = new uint8_t[msglen]();
	uint8_t* msg1 = new uint8_t[msglen]();
	uint8_t* msgb = new uint8_t[msglen]();
	for (int i = 0; i < OPRF_k; i++)
	{
		for (int j = 0; j < m / 8; j++)
		{
			uint8_t v0 = 0, v1 = 0;
			for (int l = 0; l < 8; l++)
			{
				v0 <<= 1;
				v0 |= t0[8 * j + l][i];
				v1 <<= 1;
				v1 |= t1[8 * j + l][i];
			}
			msg0[j] = v0;
			msg1[j] = v1;
		}
		for (int j = 0; j < extbits; j++)
		{
			msg0[msglen - 1] <<= 1;
			msg0[msglen - 1] |= t0[m - extbits + j][i];
			msg1[msglen - 1] <<= 1;
			msg1[msglen - 1] |= t1[m - extbits + j][i];
		}
		OT::OT(msg0, msg1, msglen, s[i], msgb);
		for (int j = 0; j < m / 8; j++)
		{
			uint8_t v = msgb[j];
			for (int l = 7; l >= 0; l--)
			{
				q[8 * j + l][i] = v & 1;
				v >>= 1;
			}
		}
		for (int j = 0; j < extbits; j++)
		{
			q[m - 1 - j][i] = msgb[msglen - 1] & 1;
			msgb[msglen - 1] >>= 1;
		}
	}
	delete[] msg0;
	delete[] msg1;
	delete[] msgb;

	outX.resize(m);
	outsY.resize(m);

	// The receiver computes outX
	for (int i = 0; i < m; i++)
		outX[i] = OPRF_H(t0[i]);

	// The sender computes outsY
	for (int i = 0; i < m; i++)
	{
		int sz = (int)elesY[i].size();
		outsY[i].resize(sz);
		for (int j = 0; j < sz; j++)
		{
			bitset<OPRF_k> bs = (PRC(i, elesY[i][j]) & s) ^ q[i];
			outsY[i][j] = OPRF_H(bs);
		}
	}
	delete[] t0;
	delete[] t1;
	delete[] q;
}
