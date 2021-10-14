#pragma once
#include<cinttypes>

class MortonCode64
{
public:
	MortonCode64() {};
	MortonCode64(int32_t x, int32_t y, int32_t z);
	MortonCode64(uint32_t x, uint32_t y, uint32_t z);
	uint64_t SplitBy3Bits21(int32_t x);
	bool operator<(const MortonCode64 rhs) const { return data < rhs.data; }
	bool operator>(const MortonCode64 rhs) const { return data > rhs.data; }
private:
	uint64_t data;
};
