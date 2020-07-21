#pragma once
#include <vector>
static const int PSI_empty = -1;

void PSI(const std::vector<uint64_t>& AliceSet, const std::vector<uint64_t>& BobSet, std::vector<int>& X_Indices, std::vector<bool>& indicator1, std::vector<bool>& indicator2);

void PSIwithPayload(const std::vector<uint64_t>& AliceSet, const std::vector<uint64_t>& BobSet, std::vector<uint16_t>& BobPayload, std::vector<int>& AliceIndicesHashed,
	std::vector<uint16_t>& paylodsOut1, std::vector<uint16_t>& paylodsOut2);

// payloadsIn1 is Alice's share
void PSIwithSharedPayload (const std::vector<uint64_t>& AliceSet, const std::vector<uint64_t>& BobSet,
	const std::vector<uint16_t>& payloadsIn1, const std::vector<uint16_t>& paylodsIn2,
	std::vector<int>& AliceIndicesHashed, std::vector<uint16_t>& payloadsOut1, std::vector<uint16_t>& payloadsOut2);