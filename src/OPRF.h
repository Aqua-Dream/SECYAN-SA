#include <cstdint>
#include <vector>

// Efficient Batched Oblivious PRF with Applications to Private Set Intersection
// this setting works only for 128-bit kappa

void batch_OPRF(std::vector<uint32_t>& eleX, std::vector<std::vector<uint32_t>>& elesY,
	std::vector<uint64_t>& outX, std::vector<std::vector<uint64_t>>& outsY);