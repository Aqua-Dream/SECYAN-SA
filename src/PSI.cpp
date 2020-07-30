#include "PSI.h"
#include "poly.h"
#include "OPRF.h"
#include "MurmurHash3.h"
#include "PRNG.h"
#include "utils.h"
#include "definitions.h"
#include "OEP.h"
#include "OT.h"
#include <random>
#include <algorithm>
#include <cstdint>

using namespace std;


// In the hash arrays, we use the first 3 elements as values of 3 hash functions that assgin elements to different bins
// We use the last element as the value of second level hashing that maps elements to domain 2^32, which is only expected to have no collision for each bin
inline void singleHash(int64_t ele, uint32_t seed, uint32_t* outarr)
{
	MurmurHash3_x64_128((char*)&ele, 8, seed, outarr);
}

void PSI::AliceCuckooHash(uint32_t** AliceHashArrs, int threshold)
{
	int totalThreshold = threshold * AliceSetSize;
	AliceIndicesHashed = new int[bucketSize];
	fill_n(AliceIndicesHashed, bucketSize, empty);
	uniform_int_distribution<> hash_id_dist(0, 2);
	for (int i = 0; i < AliceSetSize; i++)
	{
		int index = i;
		int bin_id;
		uint8_t hash_id;
		bool finished = false;
		while (1)
		{
			for (hash_id = 0; hash_id < 3; hash_id++)
			{
				bin_id = AliceHashArrs[index][hash_id] % bucketSize;
				if (AliceIndicesHashed[bin_id] == empty)
				{
					finished = true;
					AliceIndicesHashed[bin_id] = index;
					break;
				}
			}
			if (finished)
				break;
			if (totalThreshold-- < 0)
				throw "Cuckoo hash threshold exceeded!";
			hash_id = hash_id_dist(gRNG.rng);
			bin_id = AliceHashArrs[index][hash_id] % bucketSize;
			swap(index, AliceIndicesHashed[bin_id]);
		}
	}
}

void PSI::BobSimpleHash(uint32_t** BobHashArrs)
{
	int locations[3];
	BobIndexVectorsHashed = new vector<int>[bucketSize];
	for (int i = 0; i < BobSetSize; i++)
	{
		for (uint8_t hash_id = 0; hash_id < 3; hash_id++)
		{
			int loc = BobHashArrs[i][hash_id] % bucketSize;
			locations[hash_id] = loc;
			bool inserted = false;
			for (int j = 0; j < hash_id; j++)
				if (locations[j] == loc)
					inserted = true;
			if (!inserted)
				BobIndexVectorsHashed[loc].push_back(i);
		}
	}
}

PSI::PSI(const vector<uint64_t>& AliceSet, const vector<uint64_t>& BobSet)
{
	AliceSetSize = (int)AliceSet.size();
	this->AliceSet = new uint64_t[AliceSetSize];
	copy(AliceSet.begin(), AliceSet.end(), this->AliceSet);
	BobSetSize = (int)BobSet.size();
	this->BobSet = new uint64_t[BobSetSize];
	copy(BobSet.begin(), BobSet.end(), this->BobSet);
	if (AliceSetSize < 30 && BobSetSize < 30)
		throw "PSI set size too small!"; // In this case faster algorithm is required (a simple GC may be enough)
	bucketSize = max((int)(1.27 * AliceSetSize), 1 + BobSetSize / 256);
	int logBucketSize = log2_32(bucketSize);
	if (logBucketSize >= 61 - sigma)
		throw "PSI set size too large!";
	int numBinsInMegabin = max((int)(bucketSize * log2(bucketSize) / BobSetSize), 1);
	numMegabins = (bucketSize + numBinsInMegabin - 1) / numBinsInMegabin;
	int m = 3 * BobSetSize;
	int n = numMegabins;
	double logn = 0;
	if (n > 1)
		logn = log2(n);
	if (m > n * logn * 4) // Throw m balls into n bins, see  "Balls into Bins" A Simple and Tight Analysis
		megaBinLoad = (double)m / n + sqrt(2 * m / n * logn);
	else
		megaBinLoad = 1.41 * m / n + 1.04 * logn;
	if (n > 1)
		megaBinLoad += sigma / logn * (1 - 1.0 / n);
	else
		megaBinLoad = m;
	// In any case, the load is at most 3*256+20=778
	gamma = sigma + logBucketSize;
	prepare();
	computeIndicators();
}

void PSI::prepare()
{
	uint32_t seed = gRNG.nextUInt32();
	addComm(32); // Sends the seed for hash

	// Alice builds hash table
	uint32_t** AliceHashArrs = new uint32_t * [AliceSetSize];
	for (int i = 0; i < AliceSetSize; i++)
	{
		AliceHashArrs[i] = new uint32_t[4];
		singleHash(AliceSet[i], seed, AliceHashArrs[i]);
	}

	AliceCuckooHash(AliceHashArrs);
	cuckooTable.resize(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		int index = AliceIndicesHashed[i];
		cuckooTable[i] = index == empty ? empty : AliceSet[index];
	}

	// Bob builds hash table
	uint32_t** BobHashArrs = new uint32_t * [BobSetSize];
	for (int i = 0; i < BobSetSize; i++)
	{
		BobHashArrs[i] = new uint32_t[4];
		singleHash(BobSet[i], seed, BobHashArrs[i]);
	}
	BobSimpleHash(BobHashArrs);
	simpleTable.resize(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		simpleTable[i].reserve(BobIndexVectorsHashed[i].size());
		for (int index : BobIndexVectorsHashed[i])
			simpleTable[i].push_back(BobSet[index]);
	}

	// OPRF

	batch_OPRF(cuckooTable, simpleTable, encCuckooTable, encSimpleTable);
	for (int i = 0; i < AliceSetSize; i++)
		delete[] AliceHashArrs[i];
	for (int i = 0; i < BobSetSize; i++)
		delete[] BobHashArrs[i];
	delete[] AliceHashArrs;
	delete[] BobHashArrs;

}

inline uint64_t PSI_combine(uint64_t v, uint64_t j)
{
	return (v & 0x1ffffffffff) | (j << 40);
}

void PSI::computeIndicators()
{
	// polynomial communication
	uint64_t* pointX = new uint64_t[megaBinLoad];
	uint64_t* pointY = new uint64_t[megaBinLoad];
	uint64_t* coeff = new uint64_t[megaBinLoad];
	uint64_t* AliceT = new uint64_t[bucketSize]();
	uint64_t* BobT = new uint64_t[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		BobT[i] = gRNG.nextUInt64();

	int remainder = bucketSize % numMegabins;
	int division = bucketSize / numMegabins;
	int startBinId = 0;
	for (int i = 0; i < numMegabins; i++)
	{
		int endBinId = startBinId + division + (i < remainder);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX[pointId] = PSI_combine(simpleTable[j][k], j);
				pointY[pointId] = (encSimpleTable[j][k] ^ BobT[j]) & poly_modulus;
				pointId++;
			}
		}
		if (pointId > megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX[pointId] = poly_modulus - pointId;
			pointY[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX, pointY, megaBinLoad, coeff);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
			if (AliceIndicesHashed[j] != empty)
				AliceT[j] = poly_eval(coeff, PSI_combine(cuckooTable[j], j), megaBinLoad) ^ encCuckooTable[j];
		startBinId = endBinId;
	}

	indicator1 = new bool[bucketSize];
	indicator2 = new bool[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		equalityTest(AliceT[i], BobT[i], gamma, indicator1[i], indicator2[i]);
	delete[] pointX;
	delete[] pointY;
	delete[] coeff;
	delete[] AliceT;
	delete[] BobT;
}

void PSI::getIndicators(std::vector<bool>& indicator1, std::vector<bool>& indicator2)
{
	indicator1.assign(this->indicator1, this->indicator1 + bucketSize);
	indicator2.assign(this->indicator2, this->indicator2 + bucketSize);
}

void PSI::getAliceIndicesHashed(std::vector<int>& AliceIndicesHashed)
{
	AliceIndicesHashed.assign(this->AliceIndicesHashed, this->AliceIndicesHashed + bucketSize);
}

void PSI::sendPayloads(const vector<uint32_t>& BobPayloads, vector<uint32_t>& paylaodsOut1, vector<uint32_t>& payloadsOut2)
{
	vector< vector<uint32_t> > simpleValue(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		simpleTable[i].reserve(BobIndexVectorsHashed[i].size());
		for (int index : BobIndexVectorsHashed[i]) {
			simpleValue[i].push_back(BobPayloads[index]);
		}
	}
	// polynomial communication
	uint64_t* pointX2 = new uint64_t[megaBinLoad];
	uint64_t* pointY2 = new uint64_t[megaBinLoad];
	uint64_t* coeff2 = new uint64_t[megaBinLoad];
	uint64_t* AliceT2 = new uint64_t[bucketSize]();
	uint64_t* BobT2 = new uint64_t[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		BobT2[i] = gRNG.nextUInt64();

	int remainder = bucketSize % numMegabins;
	int division = bucketSize / numMegabins;
	int startBinId = 0;
	for (int i = 0; i < numMegabins; i++)
	{
		int endBinId = startBinId + division + (i < remainder);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX2[pointId] = PSI_combine(simpleTable[j][k], j);
				pointY2[pointId] = (encSimpleTable[j][k] ^ (BobPayloads[BobIndexVectorsHashed[j][k]] - BobT2[j])) & poly_modulus;
				pointId++;
			}
		}
		if (pointId > megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX2[pointId] = poly_modulus - pointId;
			pointY2[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX2, pointY2, megaBinLoad, coeff2);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
		{
			if (AliceIndicesHashed[j] != empty)
				AliceT2[j] = poly_eval(coeff2, PSI_combine(cuckooTable[j], j), megaBinLoad) ^ encCuckooTable[j];
		}
		startBinId = endBinId;
	}
	paylaodsOut1.resize(bucketSize);
	payloadsOut2.resize(bucketSize);
	for (int i = 0; i < bucketSize; ++i) {
		paylaodsOut1[i] = (uint32_t)AliceT2[i];
		payloadsOut2[i] = (uint32_t)BobT2[i];
	}
	delete[] pointX2;
	delete[] pointY2;
	delete[] coeff2;
	delete[] AliceT2;
	delete[] BobT2;
}

void PSI::sendPayloads(const vector<uint16_t>& BobPayloads, vector<uint16_t>& payloadsOut1, vector<uint16_t>& payloadsOut2, bool applyIndicator)
{
	vector<uint32_t> BobPayloads_32, payloadsOut1_32, payloadsOut2_32;
	BobPayloads_32.assign(BobPayloads.begin(), BobPayloads.end());
	sendPayloads(BobPayloads_32, payloadsOut1_32, payloadsOut2_32);
	payloadsOut1.assign(payloadsOut1_32.begin(), payloadsOut1_32.end());
	payloadsOut2.assign(payloadsOut2_32.begin(), payloadsOut2_32.end());
	if (applyIndicator)
		for (int i = 0; i < bucketSize; i++)
			oblivExtTransfer(0, 0, payloadsOut1[i], payloadsOut2[i], indicator1[i], indicator2[i], payloadsOut1[i], payloadsOut2[i]);

}

void PSI::sendSharedPayloads(const vector<uint16_t>& payloadsIn1, const vector<uint16_t>& payloadsIn2,
	vector<uint16_t>& payloadsOut1, vector<uint16_t>& payloadsOut2)
{
	int extendShareSize = BobSetSize + bucketSize;

	// each party extend the shares
	vector<uint16_t> extendValueShare1 = payloadsIn1;
	extendValueShare1.resize(extendShareSize, 0);
	vector<uint16_t> extendValueShare2 = payloadsIn2;
	extendValueShare2.resize(extendShareSize, 0);

	vector<uint32_t> rp1(extendShareSize), invrp1(extendShareSize);
	for (int i = 0; i < extendShareSize; ++i)
		rp1[i] = i;
	shuffle(rp1.begin(), rp1.end(), gRNG.rng);
	for (int i = 0; i < extendShareSize; ++i)
		invrp1[rp1[i]] = i;

	int* indices1 = new int[extendShareSize];
	uint16_t* values = new uint16_t[extendShareSize];
	uint16_t* z1 = new uint16_t[extendShareSize];
	uint16_t* z2 = new uint16_t[extendShareSize];

	for (int i = 0; i < extendShareSize; ++i) {
		values[i] = extendValueShare2[i];
		indices1[i] = rp1[i];
	}
	oblivPermutation(indices1, values, extendShareSize, z1, z2);
	for (int i = 0; i < extendShareSize; ++i)
		z1[i] += extendValueShare1[rp1[i]];

	vector<uint32_t> AliceRev, BobRev;
	sendPayloads(invrp1, AliceRev, BobRev);

	// modify PSI, if indicator1[i] ^ indicator2[i] = true, remain (AliceRev, BobRev); else change AliceRev + BobRev = 0
	for (int i = 0; i < bucketSize; ++i)
		oblivExtTransfer(0, invrp1[BobSetSize + i], AliceRev[i], BobRev[i], indicator1[i], indicator2[i], AliceRev[i], BobRev[i]);

	vector<uint32_t> ki(bucketSize);
	addComm(16 * bucketSize);
	// Reveal ki to Alice
	for (int i = 0; i < bucketSize; ++i)
		ki[i] = AliceRev[i] + BobRev[i];

	int* indices2 = new int[bucketSize];
	uint16_t* zz1 = new uint16_t[bucketSize];
	uint16_t* zz2 = new uint16_t[bucketSize];
	copy(ki.begin(), ki.end(), indices2);
	OEP(indices2, extendShareSize, bucketSize, z2, zz1, zz2);
	for (int i = 0; i < bucketSize; ++i)
		zz1[i] += z1[ki[i]];
	payloadsOut1.assign(zz1, zz1 + bucketSize);
	payloadsOut2.assign(zz2, zz2 + bucketSize);

	delete[] values;
	delete[] indices1;
	delete[] z1;
	delete[] z2;
	delete[] indices2;
	delete[] zz1;
	delete[] zz2;
}

PSI::~PSI()
{
	delete[] AliceSet;
	delete[] BobSet;
	delete[] AliceIndicesHashed;
	delete[] BobIndexVectorsHashed;
	delete[] indicator1;
	delete[] indicator2;
}