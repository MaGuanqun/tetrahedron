#include"morton.h"
#include <cassert>
uint64_t MortonCode64::SplitBy3Bits21(int32_t x)
{
    int64_t r = x & 0x1fffff;//only use first 21
    //uint64_t r = x;

    r = (r | r << 32) & 0x1f00000000ffff;//0000000000011111000000000000000000000000000000001111111111111111
    r = (r | r << 16) & 0x1f0000ff0000ff;//0000000000011111000000000000000011111111000000000000000011111111
    r = (r | r << 8) & 0x100f00f00f00f00f;//0001000000001111000000001111000000001111000000001111000000001111
    r = (r | r << 4) & 0x10c30c30c30c30c3;//0001000011000011000011000011000011000011000011000011000011000011
    r = (r | r << 2) & 0x1249249249249249;//0001001001001001001001001001001001001001001001001001001001001001
    return r;
}


MortonCode64::MortonCode64(int32_t x, int32_t y, int32_t z)
{
    assert(-(1 << 20) <= x && x < (1 << 20));
    assert(-(1 << 20) <= y && y < (1 << 20));
    assert(-(1 << 20) <= z && z < (1 << 20));
    // move sign bit to bit 20
    x = (x & 0x80000000) >> 11 | (x & 0x0fffff);
    y = (y & 0x80000000) >> 11 | (y & 0x0fffff);
    z = (z & 0x80000000) >> 11 | (z & 0x0fffff);
    data = SplitBy3Bits21(x) | SplitBy3Bits21(y) << 1 | SplitBy3Bits21(z) << 2;
    data = data ^ 0x7000000000000000;
}

MortonCode64::MortonCode64(uint32_t x, uint32_t y, uint32_t z)
{
    data = (SplitBy3Bits21(x) | SplitBy3Bits21(y) << 1 | SplitBy3Bits21(z) << 2);
    //data |= 0x7000000000000000;
}



