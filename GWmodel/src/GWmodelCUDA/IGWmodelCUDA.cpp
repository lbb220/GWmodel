#include "IGWmodelCUDA.h"
#include "CGWmodelCUDA.h"


IGWmodelCUDA * GWCUDA_Create(int N, int K, bool rp_given, int n, bool dm_given)
{
	return new CGWmodelCUDA(N, K, rp_given, n, dm_given);
}

void GWCUDA_Del(IGWmodelCUDA * pInstance)
{
	delete pInstance;
}

bool gwr_reg_cuda(IGWmodelCUDA * pInstance, bool hatmatrix, double p, double theta, bool longlat, double bw, int kernel, bool adaptive, int groupl, int gpuID)
{
	return pInstance->Regression(hatmatrix, p, theta, longlat, bw, kernel, adaptive, groupl, gpuID);
}
