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


int bucketSize, AliceSetSize, BobSetSize, numBinsInMegabin, numMegabins, megaBinLoad, PSI_gamma;
uint32_t** AliceHashArrs;
uint32_t** BobHashArrs;

void updateSizes()
{
	bucketSize = max((int)(1.27 * AliceSetSize), 1 + BobSetSize / 1024);
	int logBucketSize = log2_32(bucketSize);
	if (logBucketSize >= 20)
		throw "PSI set size too large!";
	numBinsInMegabin = max(bucketSize * logBucketSize / BobSetSize, 1);

	numMegabins = (bucketSize + numBinsInMegabin - 1) / numBinsInMegabin;
	int m = 3 * BobSetSize;
	int n = numMegabins - 1; // The load of the last mega bin is not ballanced
	if (m > n * log2(n + 1) * 4) // Throw m balls into n bins, see  "Balls into Bins" ¡ª A Simple and Tight Analysis
		megaBinLoad = (double)m / n + sqrt(2 * m / n * log2(n));
	else
		megaBinLoad = 1.41 * m / n + 1.04 * log2(n);
	megaBinLoad += sigma / log2(sigma);
	PSI_gamma = min(61, (sigma + logBucketSize));
}

// In the hash arrays, we use the first 3 elements as values of 3 hash functions that assgin elements to different bins
// We use the last element as the value of second level hashing that maps elements to domain 2^32, which is only expected to have no collision for each bin
inline void singleHash(int64_t ele, uint32_t seed, uint32_t* outarr)
{
	MurmurHash3_x64_128((char*)&ele, 8, seed, outarr);
}

void AliceCuckooHash(vector<int>& AliceIndicesHashed, int threshold = 3)
{
	int totalThreshold = threshold * AliceSetSize;
	AliceIndicesHashed.assign(bucketSize, PSI_empty);
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
				if (AliceIndicesHashed[bin_id] == PSI_empty)
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

void BobSimpleHash(vector< vector<int> >& BobIndexVectorsHashed)
{
	int locations[3];
	BobIndexVectorsHashed.resize(bucketSize);
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

void PSI(const vector<uint64_t>& AliceSet, const vector<uint64_t>& BobSet, vector<int>& AliceIndicesHashed, vector<bool>& indicator1, vector<bool>& indicator2)
{
	AliceSetSize = (int)AliceSet.size();
	BobSetSize = (int)BobSet.size();
	updateSizes();

	uint32_t seed = gRNG.nextUInt32();
	addComm(32); // Sends the seed for hash

	// Alice builds hash table
	AliceHashArrs = new uint32_t * [AliceSetSize];
	for (int i = 0; i < AliceSetSize; i++)
	{
		AliceHashArrs[i] = new uint32_t[4];
		singleHash(AliceSet[i], seed, AliceHashArrs[i]);
	}

	AliceCuckooHash(AliceIndicesHashed);
	vector<uint32_t> cuckooTable(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		int index = AliceIndicesHashed[i];
		if (index == PSI_empty)
			cuckooTable[i] = PSI_empty;
		else
			cuckooTable[i] = AliceHashArrs[index][3];
	}

	// Bob builds hash table
	BobHashArrs = new uint32_t * [BobSetSize];
	for (int i = 0; i < BobSetSize; i++)
	{
		BobHashArrs[i] = new uint32_t[4];
		singleHash(BobSet[i], seed, BobHashArrs[i]);
	}
	vector<vector<int>> BobIndexVectorsHashed;
	BobSimpleHash(BobIndexVectorsHashed);
	vector< vector<uint32_t> > simpleTable(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		simpleTable[i].reserve(BobIndexVectorsHashed[i].size());
		for (int index : BobIndexVectorsHashed[i])
			simpleTable[i].push_back(BobHashArrs[index][3]);
	}

	// OPRF
	vector<uint64_t> encCuckooTable;
	vector< vector<uint64_t>> encSimpleTable;
	batch_OPRF(cuckooTable, simpleTable, encCuckooTable, encSimpleTable);

	// polynomial communication
	uint64_t* pointX = new uint64_t[megaBinLoad];
	uint64_t* pointY = new uint64_t[megaBinLoad];
	uint64_t* coeff = new uint64_t[megaBinLoad];
	uint64_t* AliceT = new uint64_t[bucketSize]();
	uint64_t* BobT = new uint64_t[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		BobT[i] = gRNG.nextUInt64();

	for (int i = 0; i < numMegabins; i++)
	{
		int startBinId = i * numBinsInMegabin;
		int endBinId = min(startBinId + numBinsInMegabin, bucketSize);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX[pointId] = simpleTable[j][k] | (((uint64_t)j) << 32);
				pointY[pointId] = (encSimpleTable[j][k] ^ BobT[j]) & poly_modulus;
				pointId++;
			}
		}
		if (pointId >= megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX[pointId] = (uint64_t)(-pointId); // Since non-dummy values has X>=0, so there is no conflict
			pointY[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX, pointY, megaBinLoad, coeff);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
		{
			if (AliceIndicesHashed[j] != PSI_empty)
				AliceT[j] = poly_eval(coeff, cuckooTable[j] | (((uint64_t)j) << 32), megaBinLoad) ^ encCuckooTable[j];
		}
	}

	indicator1.resize(bucketSize);
	indicator2.resize(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		bool z1 = 0, z2 = 0;
		equalityTest(AliceT[i], BobT[i], PSI_gamma, z1, z2);
		indicator1[i] = z1;
		indicator2[i] = z2;
	}

	for (int i = 0; i < AliceSetSize; i++)
		delete[] AliceHashArrs[i];
	for (int i = 0; i < BobSetSize; i++)
		delete[] BobHashArrs[i];
	delete[] AliceHashArrs;
	delete[] BobHashArrs;
	delete[] pointX;
	delete[] pointY;
	delete[] coeff;
	delete[] AliceT;
	delete[] BobT;
}

void PSIwithPayload(const vector<uint64_t>& AliceSet, const vector<uint64_t>& BobSet, vector<uint32_t>& BobPayload, vector<int>& AliceIndicesHashed,
	vector<uint32_t>& paylodsOut1, vector<uint32_t>& paylodsOut2,
	vector<bool>& indicator1, vector<bool>& indicator2)
{
	AliceSetSize = (int)AliceSet.size();
	BobSetSize = (int)BobSet.size();
	updateSizes();
	uint32_t seed = gRNG.nextUInt32();
	addComm(32); // Sends the seed for hash

	// Alice builds hash table
	AliceHashArrs = new uint32_t * [AliceSetSize];
	for (int i = 0; i < AliceSetSize; i++)
	{
		AliceHashArrs[i] = new uint32_t[4];
		singleHash(AliceSet[i], seed, AliceHashArrs[i]);
	}

	AliceCuckooHash(AliceIndicesHashed);
	vector<uint32_t> cuckooTable(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		int index = AliceIndicesHashed[i];
		if (index == PSI_empty)
			cuckooTable[i] = PSI_empty;
		else
			cuckooTable[i] = AliceHashArrs[index][3];
	}

	// Bob builds hash table
	BobHashArrs = new uint32_t * [BobSetSize];
	for (int i = 0; i < BobSetSize; i++)
	{
		BobHashArrs[i] = new uint32_t[4];
		singleHash(BobSet[i], seed, BobHashArrs[i]);
	}
	vector<vector<int>> BobIndexVectorsHashed;
	BobSimpleHash(BobIndexVectorsHashed);
	vector< vector<uint32_t> > simpleTable(bucketSize);
	vector< vector<uint32_t> > simpleValue(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		simpleTable[i].reserve(BobIndexVectorsHashed[i].size());
		for (int index : BobIndexVectorsHashed[i]) {
			simpleTable[i].push_back(BobHashArrs[index][3]);
			simpleValue[i].push_back(BobPayload[index]);
		}
	}

	// OPRF
	vector<uint64_t> encCuckooTable;
	vector< vector<uint64_t>> encSimpleTable;
	batch_OPRF(cuckooTable, simpleTable, encCuckooTable, encSimpleTable);

	// polynomial communication
	uint64_t* pointX1 = new uint64_t[megaBinLoad];
	uint64_t* pointY1 = new uint64_t[megaBinLoad];
	uint64_t* coeff1 = new uint64_t[megaBinLoad];
	uint64_t* AliceT1 = new uint64_t[bucketSize]();
	uint64_t* BobT1 = new uint64_t[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		BobT1[i] = gRNG.nextUInt64();

	for (int i = 0; i < numMegabins; i++)
	{
		int startBinId = i * numBinsInMegabin;
		int endBinId = min(startBinId + numBinsInMegabin, bucketSize);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX1[pointId] = simpleTable[j][k] | (((uint64_t)j) << 32);
				pointY1[pointId] = (encSimpleTable[j][k] ^ BobT1[j]) & poly_modulus;
				pointId++;
			}
		}
		if (pointId >= megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX1[pointId] = (uint64_t)(-pointId); // Since non-dummy values has X>=0, so there is no conflict
			pointY1[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX1, pointY1, megaBinLoad, coeff1);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
		{
			if (AliceIndicesHashed[j] != PSI_empty)
				AliceT1[j] = poly_eval(coeff1, cuckooTable[j] | (((uint64_t)j) << 32), megaBinLoad) ^ encCuckooTable[j];
		}
	}
	indicator1.resize(bucketSize);
	indicator2.resize(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		bool z1 = 0, z2 = 0;
		equalityTest(AliceT1[i], BobT1[i], PSI_gamma, z1, z2);
		indicator1[i] = z1;
		indicator2[i] = z2;
	}

	// polynomial communication
	uint64_t* pointX2 = new uint64_t[megaBinLoad];
	uint64_t* pointY2 = new uint64_t[megaBinLoad];
	uint64_t* coeff2 = new uint64_t[megaBinLoad];
	uint64_t* AliceT2 = new uint64_t[bucketSize]();
	uint64_t* BobT2 = new uint64_t[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		BobT2[i] = gRNG.nextUInt64();

	for (int i = 0; i < numMegabins; i++)
	{
		int startBinId = i * numBinsInMegabin;
		int endBinId = min(startBinId + numBinsInMegabin, bucketSize);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX2[pointId] = simpleTable[j][k] | (((uint64_t)j) << 32);
				pointY2[pointId] = (encSimpleTable[j][k] ^ ((uint64_t)simpleValue[j][k] - BobT2[j])) & poly_modulus;
				pointId++;
			}
		}
		if (pointId >= megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX2[pointId] = (uint64_t)(-pointId); // Since non-dummy values has X>=0, so there is no conflict
			pointY2[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX2, pointY2, megaBinLoad, coeff2);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
		{
			if (AliceIndicesHashed[j] != PSI_empty)
				AliceT2[j] = poly_eval(coeff2, cuckooTable[j] | (((uint64_t)j) << 32), megaBinLoad) ^ encCuckooTable[j];
		}
	}
	paylodsOut1.resize(bucketSize);
	paylodsOut2.resize(bucketSize);
	for (int i = 0; i < bucketSize; ++i) {
		paylodsOut1[i] = (uint32_t)AliceT2[i];
		paylodsOut2[i] = (uint32_t)BobT2[i];
	}
	for (int i = 0; i < AliceSetSize; i++)
		delete[] AliceHashArrs[i];
	for (int i = 0; i < BobSetSize; i++)
		delete[] BobHashArrs[i];
	delete[] AliceHashArrs;
	delete[] BobHashArrs;
	delete[] pointX1;
	delete[] pointY1;
	delete[] coeff1;
	delete[] AliceT1;
	delete[] BobT1;
	delete[] pointX2;
	delete[] pointY2;
	delete[] coeff2;
	delete[] AliceT2;
	delete[] BobT2;

}

void PSIwithSharedPayload(const vector<uint64_t>& AliceSet, const vector<uint64_t>& BobSet,
	const vector<uint16_t>& payloadsIn1, const vector<uint16_t>& payloadsIn2,
	vector<int>& AliceIndicesHashed, vector<uint16_t>& payloadsOut1, vector<uint16_t>& payloadsOut2)
{
	AliceSetSize = (int)AliceSet.size();
	BobSetSize = (int)BobSet.size();
	updateSizes();
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
	vector<bool> indicator1, indicator2;
	PSIwithPayload(AliceSet, BobSet, invrp1, AliceIndicesHashed, AliceRev, BobRev, indicator1, indicator2);
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


void PSIwithPayload(const vector<uint64_t>& AliceSet, const vector<uint64_t>& BobSet, vector<uint16_t>& BobPayload, vector<int>& AliceIndicesHashed,
	vector<uint16_t>& paylodsOut1, vector<uint16_t>& paylodsOut2)
{
	AliceSetSize = (int)AliceSet.size();
	BobSetSize = (int)BobSet.size();
	updateSizes();
	uint32_t seed = gRNG.nextUInt32();
	addComm(32); // Sends the seed for hash

	// Alice builds hash table
	AliceHashArrs = new uint32_t * [AliceSetSize];
	for (int i = 0; i < AliceSetSize; i++)
	{
		AliceHashArrs[i] = new uint32_t[4];
		singleHash(AliceSet[i], seed, AliceHashArrs[i]);
	}

	AliceCuckooHash(AliceIndicesHashed);
	vector<uint32_t> cuckooTable(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		int index = AliceIndicesHashed[i];
		if (index == PSI_empty)
			cuckooTable[i] = PSI_empty;
		else
			cuckooTable[i] = AliceHashArrs[index][3];
	}

	// Bob builds hash table
	BobHashArrs = new uint32_t * [BobSetSize];
	for (int i = 0; i < BobSetSize; i++)
	{
		BobHashArrs[i] = new uint32_t[4];
		singleHash(BobSet[i], seed, BobHashArrs[i]);
	}
	vector<vector<int>> BobIndexVectorsHashed;
	BobSimpleHash(BobIndexVectorsHashed);
	vector< vector<uint32_t> > simpleTable(bucketSize);
	vector< vector<uint32_t> > simpleValue(bucketSize);
	for (int i = 0; i < bucketSize; i++)
	{
		simpleTable[i].reserve(BobIndexVectorsHashed[i].size());
		for (int index : BobIndexVectorsHashed[i]) {
			simpleTable[i].push_back(BobHashArrs[index][3]);
			simpleValue[i].push_back(BobPayload[index]);
		}
	}

	// OPRF
	vector<uint64_t> encCuckooTable;
	vector< vector<uint64_t>> encSimpleTable;
	batch_OPRF(cuckooTable, simpleTable, encCuckooTable, encSimpleTable);

	// polynomial communication
	uint64_t* pointX1 = new uint64_t[megaBinLoad];
	uint64_t* pointY1 = new uint64_t[megaBinLoad];
	uint64_t* coeff1 = new uint64_t[megaBinLoad];
	uint64_t* AliceT1 = new uint64_t[bucketSize]();
	uint64_t* BobT1 = new uint64_t[bucketSize];
	for (int i = 0; i < bucketSize; i++)
		BobT1[i] = gRNG.nextUInt64();

	for (int i = 0; i < numMegabins; i++)
	{
		int startBinId = i * numBinsInMegabin;
		int endBinId = min(startBinId + numBinsInMegabin, bucketSize);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX1[pointId] = simpleTable[j][k] | (((uint64_t)j) << 32);
				pointY1[pointId] = (encSimpleTable[j][k] ^ BobT1[j]) & poly_modulus;
				pointId++;
			}
		}
		if (pointId >= megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX1[pointId] = (uint64_t)(-pointId); // Since non-dummy values has X>=0, so there is no conflict
			pointY1[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX1, pointY1, megaBinLoad, coeff1);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
		{
			if (AliceIndicesHashed[j] != PSI_empty)
				AliceT1[j] = poly_eval(coeff1, cuckooTable[j] | (((uint64_t)j) << 32), megaBinLoad) ^ encCuckooTable[j];
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

	for (int i = 0; i < numMegabins; i++)
	{
		int startBinId = i * numBinsInMegabin;
		int endBinId = min(startBinId + numBinsInMegabin, bucketSize);
		int pointId = 0;
		for (int j = startBinId; j < endBinId; j++)
		{
			for (int k = 0; k < simpleTable[j].size(); k++)
			{
				pointX2[pointId] = simpleTable[j][k] | (((uint64_t)j) << 32);
				pointY2[pointId] = (encSimpleTable[j][k] ^ ((uint64_t)simpleValue[j][k] - BobT2[j])) & poly_modulus;
				pointId++;
			}
		}
		if (pointId >= megaBinLoad)
			throw "Mega bin load not enough!";
		while (pointId < megaBinLoad)
		{
			pointX2[pointId] = (uint64_t)(-pointId); // Since non-dummy values has X>=0, so there is no conflict
			pointY2[pointId] = gRNG.nextUInt64() & poly_modulus;
			pointId++;
		}
		interpolate(pointX2, pointY2, megaBinLoad, coeff2);
		addComm(61 * megaBinLoad); // Sends polynomials
		for (int j = startBinId; j < endBinId; j++)
		{
			if (AliceIndicesHashed[j] != PSI_empty)
				AliceT2[j] = poly_eval(coeff2, cuckooTable[j] | (((uint64_t)j) << 32), megaBinLoad) ^ encCuckooTable[j];
		}
	}
	paylodsOut1.resize(bucketSize);
	paylodsOut2.resize(bucketSize);

	for (int i = 0; i < bucketSize; i++)
	{
		bool z1 = 0, z2 = 0;
		equalityTest(AliceT1[i], BobT1[i], PSI_gamma, z1, z2);
		oblivExtTransfer(0, 0, (uint16_t)AliceT2[i], (uint16_t)BobT2[i], z1, z2, paylodsOut1[i], paylodsOut2[i]);
	}

	for (int i = 0; i < AliceSetSize; i++)
		delete[] AliceHashArrs[i];
	for (int i = 0; i < BobSetSize; i++)
		delete[] BobHashArrs[i];
	delete[] AliceHashArrs;
	delete[] BobHashArrs;
	delete[] pointX1;
	delete[] pointY1;
	delete[] coeff1;
	delete[] AliceT1;
	delete[] BobT1;
	delete[] pointX2;
	delete[] pointY2;
	delete[] coeff2;
	delete[] AliceT2;
	delete[] BobT2;

}