#pragma once
#include <armadillo>

#include "IGWmodelCUDA.h"

using namespace arma;

class CGWmodelCUDA : public IGWmodelCUDA
{
private:
	mat x;
	vec y;
	mat dp;
	mat rp;
	mat dMat;
	bool rp_given;
	bool dm_given;
	mat betas;
	mat betasSE;
	vec s_hat;
	vec qdiag;
public:
	CGWmodelCUDA();
	CGWmodelCUDA(int N, int K, bool rp_given, int n, bool dm_given);
	~CGWmodelCUDA();

	virtual void SetX(int i, int k, double value);
	virtual void SetY(int i, double value);
	virtual void SetDp(int i, double u, double v);
	virtual void SetRp(int i, double u, double v);
	virtual void SetDmat(int i, int j, double value);

	virtual double GetBetas(int i, int k);
	virtual double GetBetasSE(int i, int k);
	virtual double GetShat1();
	virtual double GetShat2();
	virtual double GetQdiag(int i);


	virtual bool Regression(
		bool hatmatrix,
		double p, double theta, bool longlat,
		double bw, int kernel, bool adaptive,
		int groupl, int gpuID
	);

	virtual double CV(
		double p, double theta, bool longlat,
		double bw, int kernel, bool adaptive,
		int groupl, int gpuID
	);

};
