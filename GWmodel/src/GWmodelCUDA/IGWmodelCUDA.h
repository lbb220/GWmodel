#pragma once

#ifdef WIN32
#ifdef CREATDLL_EXPORTS
#define GWMODELCUDA_API __declspec(dllexport)
#else
#define GWMODELCUDA_API __declspec(dllimport)
#endif // CREATDLL_EXPORTS
#else
#define GWMODELCUDA_API 
#endif

class GWMODELCUDA_API IGWmodelCUDA 
{
public:
	virtual void SetX(int i, int k, double value) = 0;
	virtual void SetY(int i, double value) = 0;
	virtual void SetDp(int i, double u, double v) = 0;
	virtual void SetRp(int i, double u, double v) = 0;
	virtual void SetDmat(int i, int j, double value) = 0;

	virtual double GetBetas(int i, int k) = 0;
	virtual double GetBetasSE(int i, int k) = 0;
	virtual double GetShat1() = 0;
	virtual double GetShat2() = 0;
	virtual double GetQdiag(int i) = 0;


	virtual bool Regression(
		bool hatmatrix,
		double p, double theta, bool longlat,
		double bw, int kernel, bool adaptive,
		int groupl, int gpuID
	) = 0;

	virtual double CV(
		double p, double theta, bool longlat,
		double bw, int kernel, bool adaptive,
		int groupl, int gpuID
	) = 0;
};

extern "C" GWMODELCUDA_API IGWmodelCUDA* GWCUDA_Create(int N, int K, bool rp_given, int n, bool dm_given);
extern "C" GWMODELCUDA_API void GWCUDA_Del(IGWmodelCUDA* pInstance);
