#pragma once
namespace Basic {
	template<typename T>
	void cumulativeSum(T* cell_num, T size, T* prefix_sum)
	{
		prefix_sum[0] = 0;
		for (T i = 0; i < size; ++i) {			
			prefix_sum[i + 1] = prefix_sum[i] + cell_num[i];
		}
	}
	template<typename T>
	void cumulativeSumDoubleRecord(T* cell_num, T size, T* prefix_sum)	//cell num store double actual num
	{
		prefix_sum[0] = 0;
		for (T i = 0; i < size; ++i) {
			prefix_sum[i + 1] = prefix_sum[i] + (cell_num[i]>>1);
		}
	}
};
