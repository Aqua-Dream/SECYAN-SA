#pragma once
#include <vector>

class PSI
{
public:
	const int empty = -1;
	PSI(const std::vector<uint64_t>& AliceSet, const std::vector<uint64_t>& BobSet);
	~PSI();
	void sendPayloads(const std::vector<uint16_t>& BobPayloads, std::vector<uint16_t>& payloadsOut1, std::vector<uint16_t>& payloadsOut2, bool applyIndicator = true);
	// payloadsIn1 is Alice's share
	void sendSharedPayloads(const std::vector<uint16_t>& payloadsIn1, const std::vector<uint16_t>& paylodsIn2,
		std::vector<uint16_t>& payloadsOut1, std::vector<uint16_t>& payloadsOut2);
	void getIndicators(std::vector<bool>& indicator1, std::vector<bool>& indicator2);
	void getAliceIndicesHashed(std::vector<int> &AliceIndicesHashed);
private:
	int bucketSize, AliceSetSize, BobSetSize, numMegabins, megaBinLoad, gamma;
	uint64_t* AliceSet, * BobSet;
	std::vector<uint64_t> cuckooTable, encCuckooTable;
	std::vector< std::vector<uint64_t> > simpleTable, encSimpleTable;
	int* AliceIndicesHashed;
	std::vector<int>* BobIndexVectorsHashed;
	bool* indicator1;
	bool* indicator2;
	void AliceCuckooHash(uint32_t** AliceHashArrs, int threshold = 3);
	void BobSimpleHash(uint32_t** BobHashArrs);
	void prepare();
	void computeIndicators();
	void sendPayloads(const std::vector<uint32_t>& BobPayloads, std::vector<uint32_t>& payloadsOut1, std::vector<uint32_t>& payloadsOut2);

};