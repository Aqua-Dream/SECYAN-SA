#include "OT.h"
#include "blake3.h"
#include "utils.h"
#include "definitions.h"
#include "PRNG.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

using namespace std;



namespace OT
{
	class bitset
	{
	public:
		bitset operator^(const bitset& bs)
		{
			bitset ret;
			for (int i = 0; i < bufsize; i++)
				ret.buffer[i] = this->buffer[i] ^ bs.buffer[i];
			return ret;
		}
		bitset operator&(const bitset& bs)
		{
			bitset ret;
			for (int i = 0; i < bufsize; i++)
				ret.buffer[i] = this->buffer[i] & bs.buffer[i];
			return ret;
		}
		bool operator[](int k)
		{
			return buffer[k / 64] & (1ull << ((~k) & 63));
		}
		uint32_t* uInt32Pt()
		{
			return (uint32_t*)buffer;
		}
		static bitset random()
		{
			bitset ret;
			for (int i = 0; i < bufsize; i++)
				ret.buffer[i] = gRNG.nextUInt64();
			return ret;
		}
		static bitset fromUInt64Arr(uint64_t* buffer)
		{
			bitset bs;
			copy(buffer, buffer + bufsize, bs.buffer);
			return bs;
		}
		static bitset fromZZ(const NTL::ZZ& zz)
		{
			bitset ret;
			NTL::BytesFromZZ((uint8_t*)ret.buffer, zz, kappa / 8);
			return ret;
		}
		static bitset fromZZp(const NTL::ZZ_p& zzp)
		{
			return fromZZ(NTL::rep(zzp));
		}
		NTL::ZZ toZZ()
		{
			return NTL::ZZFromBytes((uint8_t*)buffer, kappa / 8);
		}
		bitset randomOracle() // for base OT
		{
			bitset ret;
			blake3_hasher hasher;
			blake3_hasher_init(&hasher);
			blake3_hasher_update(&hasher, buffer, kappa / 8);
			blake3_hasher_finalize(&hasher, (uint8_t*)ret.buffer, kappa / 8);
			return ret;
		}
		void randomOracle(uint64_t seed, uint8_t* out, int outlen) // for OT
		{
			uint64_t newbuf[bufsize + 1];
			copy(buffer, buffer + bufsize, newbuf);
			newbuf[bufsize] = seed;
			blake3_hasher hasher;
			blake3_hasher_init(&hasher);
			blake3_hasher_update(&hasher, newbuf, kappa / 8 + 8);
			blake3_hasher_finalize(&hasher, out, outlen);
		}
		bool randomOracle(uint64_t seed) // for GMW OT
		{
			uint8_t out;
			uint64_t newbuf[bufsize + 1];
			copy(buffer, buffer + bufsize, newbuf);
			newbuf[bufsize] = seed;
			blake3_hasher hasher;
			blake3_hasher_init(&hasher);
			blake3_hasher_update(&hasher, newbuf, kappa / 8 + 8);
			blake3_hasher_finalize(&hasher, &out, 1);
			return out & 1;
		}
	private:
		static const int bufsize = kappa / 64;
		uint64_t buffer[bufsize];
	};

	static bool init;
	static uint64_t index;
	NTL::ZZ P;
	NTL::ZZ_p g;
	bitset s;
	PRNG gk0[kappa], gk1[kappa], gks[kappa];

	// The messages must be 128-bit numbers
	bitset baseOT_128(bitset msg0, bitset msg1, bool b)
	{
		using namespace NTL;
		// g is the generator and P is a big prime, both public to Sender and Receiver

		// Sender randomly generates C, and send C to Receiver
		ZZ_p C = random_ZZ_p();
		// Sends C
		addComm(kappa);

		// Receiver generates random number k and computes PK0 and PK1, and send PK0 to Sender
		ZZ_p receiverPK[2];
		ZZ k = RandomBits_ZZ(kappa - 1);
		power(receiverPK[b], g, k);
		receiverPK[!b] = C / receiverPK[b];
		// Sends receiverPK[0]
		addComm(kappa);

		ZZ_p senderPK[2];
		senderPK[0] = receiverPK[0];
		senderPK[1] = C / senderPK[0];

		ZZ r[2];
		RandomBits(r[0], kappa - 1);
		RandomBits(r[1], kappa - 1);

		bitset M[2];
		M[0] = msg0;
		M[1] = msg1;

		bitset buffer[4]; // message to send
		for (int i = 0; i < 2; ++i)
		{
			buffer[2 * i] = bitset::fromZZp(power(g, r[i]));
			buffer[2 * i + 1] = bitset::fromZZp(power(senderPK[i], r[i])).randomOracle() ^ M[i];
		}
		addComm(kappa * 4); // sending message
		return bitset::fromZZ(PowerMod(buffer[2 * b].toZZ(), k, P)).randomOracle() ^ buffer[2 * b + 1];
	}

	void OT_Pre()
	{
		NTL::conv(P, "95547308683929808626692633782906550161");
		NTL::ZZ_p::init(P);
		NTL::conv(g, 65537);
		s = bitset::random();

		// Receiver randomly generates kappa pairs of kappa bits seed
		for (int i = 0; i < kappa; i++)
		{
			bitset sdk0 = bitset::random();
			bitset sdk1 = bitset::random();
			bitset sdks = baseOT_128(sdk0, sdk1, s[i]);
			seed_seq seed0(sdk0.uInt32Pt(), sdk0.uInt32Pt() + kappa / 32);
			gk0[i].seed(seed0);
			seed_seq seed1(sdk1.uInt32Pt(), sdk1.uInt32Pt() + kappa / 32);
			gk1[i].seed(seed1);
			seed_seq seeds(sdks.uInt32Pt(), sdks.uInt32Pt() + kappa / 32);
			gks[i].seed(seeds);
		}
	}

	void OT_Common(bool b, bitset& t, bitset& q, bitset& qs)
	{
		if (!init)
		{
			OT_Pre();
			init = true;
		}
		if (index++ == UINT64_MAX)
			throw "OT index used up!";

		uint64_t t_arr[kappa / 64], u_arr[kappa / 64];

		for (int i = 0; i < kappa / 64; i++)
		{
			uint64_t t_arri = 0, u_arri = 0;
			for (int j = 0; j < 64; j++)
			{
				bool ti = gk0[64 * i + j].nextBit();
				t_arri = (t_arri << 1) | ti;
				bool ui = gk1[64 * i + j].nextBit() ^ ti;
				u_arri = (u_arri << 1) | ui;
			}
			if (b)
				u_arri = ~u_arri;
			t_arr[i] = t_arri;
			u_arr[i] = u_arri;
		}
		t = bitset::fromUInt64Arr(t_arr);
		bitset u = bitset::fromUInt64Arr(u_arr);

		// Sends u
		addComm(kappa);

		uint64_t gksi_arr[kappa / 64];
		for (int i = 0; i < kappa / 64; i++)
		{
			uint64_t gksi_arri = 0;
			for (int j = 0; j < 64; j++)
				gksi_arri = (gksi_arri << 1) | gks[64 * i + j].nextBit();
			gksi_arr[i] = gksi_arri;
		}

		q = (s & u) ^ bitset::fromUInt64Arr(gksi_arr);
		qs = q ^ s;
	}

	void OT(uint8_t* msg0, uint8_t* msg1, int length, bool b, uint8_t* msgb)
	{

		bitset t, q, qs;
		OT_Common(b, t, q, qs);
		uint8_t* y0 = new uint8_t[length];
		uint8_t* y1 = new uint8_t[length];
		q.randomOracle(index, y0, length);
		qs.randomOracle(index, y1, length);
		for (int i = 0; i < length; i++)
		{
			y0[i] ^= msg0[i];
			y1[i] ^= msg1[i];
		}
		// Sends buffer
		addComm(8 * length * 2);

		uint8_t* yb = b ? y1 : y0;
		t.randomOracle(index, msgb, length);
		for (int i = 0; i < length; i++)
			msgb[i] ^= yb[i];
		delete[] y0;
		delete[] y1;
	}


	// The special OT only for GMW
	void GMW_OT(bool b, bool& msg0, bool& msg1, bool& msgb)
	{
		bitset t, q, qs;
		OT_Common(b, t, q, qs);
		msg0 = q.randomOracle(index);
		msg1 = qs.randomOracle(index);
		msgb = t.randomOracle(index);
	}


	// msg0 is random and msg1==msg0^delta
	void correlatedOT(bool b, uint8_t* delta, int length, uint8_t* msg0, uint8_t* msgb)
	{
		bitset t, q, qs;
		OT_Common(b, t, q, qs);
		q.randomOracle(index, msg0, length);
		uint8_t* qs_hasharr = new uint8_t[length];
		qs.randomOracle(index, qs_hasharr, length);
		uint8_t* buffer = new uint8_t[length];
		for (int i = 0; i < length; i++)
			buffer[i] = delta[i] ^ msg0[i] ^ qs_hasharr[i];
		// Sends buffer
		addComm(8 * length);

		t.randomOracle(index, msgb, length);
		if (b)
			for (int i = 0; i < length; i++)
				msgb[i] ^= buffer[i];
		delete[] qs_hasharr;
		delete[] buffer;
	}


	void shareMulCOT(bool b, uint16_t delta, uint16_t& msg0, uint16_t& msgb)
	{
		bitset t, q, qs;
		uint16_t buffer;
		OT_Common(b, t, q, qs);
		q.randomOracle(index, (uint8_t*)&msg0, 2);
		qs.randomOracle(index, (uint8_t*)&buffer, 2);
		buffer += delta - msg0;
		addComm(16); // send the buffer

		uint16_t msgcal;
		t.randomOracle(index, (uint8_t*)&msgcal, 2);
		msgb = (b ? buffer : 0) - msgcal;
	}


}