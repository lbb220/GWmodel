#include <armadillo>

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <stdio.h>

#include <time.h>

#include "helper.h"
#include "CGWmodelCUDA.h"
#include "GWmodelKernel.h"
#ifdef ENABLE_PREVIEW
#include "preview.h"
#endif // ENABLE_PREVIEW

using namespace arma;

#define RESERVE_GPU_MEM ((long long int)2*1024*1024*1024)

CGWmodelCUDA::CGWmodelCUDA()
	: x(0, 0)
	, y((uword)0)
	, dp(0, 0)
	, rp(0, 0)
	, dMat(0, 0)
	, betas(0, 0)
	, betasSE(0, 0)
	, s_hat((uword)0)
	, qdiag((uword)0)
{
	rp_given = false;
	dm_given = false;
}

CGWmodelCUDA::CGWmodelCUDA(int N, int K, bool rp_given, int n, bool dm_given)
	: x(K, N, fill::zeros)
	, y(N, fill::zeros)
	, dp(N, 2, fill::zeros)
	, betas(0, 0)
	, betasSE(0, 0)
	, s_hat((uword)0)
	, qdiag((uword)0)
{
	this->rp_given = rp_given;
	this->dm_given = dm_given;
	if (rp_given) rp = mat(n, 2, fill::zeros);
	else rp = mat(N, 2);
	if (dm_given) dMat = mat(N, n, fill::zeros);
	else dMat = mat(0, 0);
}

CGWmodelCUDA::~CGWmodelCUDA()
{
	x.clear();
	y.clear();
	dp.clear();
	rp.clear();
	dMat.clear(); 
	rp_given = false;
	dm_given = false;
}

void CGWmodelCUDA::SetX(int i, int k, double value)
{
	x(i, k) = value;
}

void CGWmodelCUDA::SetY(int i, double value)
{
	y(i) = value;
}

void CGWmodelCUDA::SetDp(int i, double u, double v)
{
	dp(i, 0) = u;
	dp(i, 1) = v;
}

void CGWmodelCUDA::SetRp(int i, double u, double v)
{
	rp(i, 0) = u;
	rp(i, 1) = v;
}

void CGWmodelCUDA::SetDmat(int i, int j, double value)
{
	dMat(i, j) = value;
}

double CGWmodelCUDA::GetBetas(int i, int k)
{
	return betas(k, i);
}

double CGWmodelCUDA::GetBetasSE(int i, int k)
{
	return betasSE(k, i);
}

double CGWmodelCUDA::GetShat1()
{
	return s_hat(0);
}

double CGWmodelCUDA::GetShat2()
{
	return s_hat(1);
}

double CGWmodelCUDA::GetQdiag(int i)
{
	return qdiag(i);
}

bool CGWmodelCUDA::Regression(bool hatmatrix,
	double p, double theta, bool longlat,
	double bw, int kernel, bool adaptive,
	int groupl, int gpuID)
{
	if (rp_given) hatmatrix = false;
	int n = rp.n_rows, N = dp.n_rows, K = x.n_rows;
	// ====================
	// Check GPU properties
	// ====================
	int gpuCount = 0;
	checkRegCudaErrors(cudaGetDeviceCount(&gpuCount));
	if (gpuID >= gpuCount)
	{
		return false;
	}
	cudaDeviceProp devProp;
	checkRegCudaErrors(cudaGetDeviceProperties(&devProp, gpuID));
	int smNum = devProp.multiProcessorCount;
	if (groupl <= 0)
	{
		groupl = (devProp.totalGlobalMem - RESERVE_GPU_MEM) / N / K / sizeof(double) / smNum * smNum;
	}
	int maxThreads = devProp.maxThreadsPerBlock;

#ifdef PRINT_CLOCKS
	printf("GPU Device: %s. \t Group length: %d\n", devProp.name, groupl);
#endif // PRINT_CLOCKS


	// ============
	// Cuda prepare
	// ============
	cublasHandle_t handle;
	cublasCreate(&handle);

	// ====================
	// Armadillo Precompute
	// ====================
	betas = mat(K, n, fill::zeros);

	// ============
	// dp rp prepare
	// ============
	double *d_dp, *d_rp;
	checkRegCudaErrors(cudaMalloc((void **)&d_dp, sizeof(double) * N * 2));
	checkRegCudaErrors(cudaMemcpy(d_dp, dp.mem, sizeof(double) * N * 2, cudaMemcpyHostToDevice));
	if (rp_given)
	{
		checkRegCudaErrors(cudaMalloc((void **)&d_rp, sizeof(double) * n * 2));
		checkRegCudaErrors(cudaMemcpy(d_rp, rp.mem, sizeof(double) * n * 2, cudaMemcpyHostToDevice));
	}
	else d_rp = d_dp;
	if (p != 2.0 && theta != 0.0 && !longlat)
	{
		checkRegCudaErrors(gw_coordinate_rotate_cuda(d_dp, N, theta, maxThreads));
		if (rp_given)
		{
			checkRegCudaErrors(gw_coordinate_rotate_cuda(d_rp, n, theta, maxThreads));
		}
	}

	// =========================
	// Regression with hatmatrix
	// =========================
	if (hatmatrix)
	{
		betasSE = mat(K, n, fill::zeros);
		s_hat = vec(2, fill::zeros);
		qdiag = vec(N, fill::zeros);

		// ----------------
		// Group Parameters
		// ----------------
		int groups = n / groupl + (n % groupl == 0 ? 0 : 1);

		// -----------------
		// CUDA memory alloc
		// -----------------
		// Single matrix
		double *d_x, *d_y, *d_dists, *d_weights;
		checkRegCudaErrors(cudaMalloc((void **)&d_x, sizeof(double) * K * N));
		checkRegCudaErrors(cudaMalloc((void **)&d_y, sizeof(double) * N * 1));
		checkRegCudaErrors(cudaMalloc((void **)&d_dists, sizeof(double) * N * 1));
		checkRegCudaErrors(cudaMalloc((void **)&d_weights, sizeof(double) * N * 1));
		checkRegCudaErrors(cudaMemcpy(d_x, x.mem, sizeof(double) * N * K, cudaMemcpyHostToDevice));
		checkRegCudaErrors(cudaMemcpy(d_y, y.mem, sizeof(double) * N * 1, cudaMemcpyHostToDevice));
		// Matrix array
		double **p_xA = new double*[groupl], **p_yA = new double*[groupl];
		for (size_t i = 0; i < groupl; i++)
		{
			p_xA[i] = d_x;
			p_yA[i] = d_y;
		}
		double **p_xiA = new double*[N], **d_xiA;
		for (size_t i = 0; i < N; i++)
		{
			p_xiA[i] = d_x + i * K;
		}
		checkRegCudaErrors(cudaMalloc((void **)&d_xiA, sizeof(double*) * N));
		checkRegCudaErrors(cudaMemcpy(d_xiA, p_xiA, sizeof(double*) * N, cudaMemcpyHostToDevice));
		double **p_xtwA, **p_xtwxA, **p_xtwyA, **p_xtwxRA, **p_betaA;
		p_xtwA = new double*[groupl];
		p_xtwxA = new double*[groupl];
		p_xtwyA = new double*[groupl];
		p_xtwxRA = new double*[groupl];
		p_betaA = new double*[groupl];
		for (size_t i = 0; i < groupl; i++)
		{
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwA[i], sizeof(double) * K * N));
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwxA[i], sizeof(double) * K * K));
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwyA[i], sizeof(double) * K * 1));
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwxRA[i], sizeof(double) * K * K));
			checkRegCudaErrors(cudaMalloc((void **)&p_betaA[i], sizeof(double) * K * 1));
		}
		double **d_xtwxA, **d_xtwxRA;
		checkRegCudaErrors(cudaMalloc((void **)&d_xtwxA, sizeof(double*) * groupl));
		checkRegCudaErrors(cudaMalloc((void **)&d_xtwxRA, sizeof(double*) * groupl));
		checkRegCudaErrors(cudaMemcpy(d_xtwxA, p_xtwxA, sizeof(double*) * groupl, cudaMemcpyHostToDevice));
		checkRegCudaErrors(cudaMemcpy(d_xtwxRA, p_xtwxRA, sizeof(double*) * groupl, cudaMemcpyHostToDevice));
		// Info vector
		int *d_info;
		checkRegCudaErrors(cudaMalloc((void **)&d_info, sizeof(int) * n));
		// qdiag
		double *d_qdiag, *d_S, *d_C;
		checkRegCudaErrors(cudaMalloc((void **)&d_S, sizeof(double) * 1 * N));
		checkRegCudaErrors(cudaMalloc((void **)&d_C, sizeof(double) * N * K));
		checkRegCudaErrors(cudaMalloc((void **)&d_qdiag, sizeof(double) * N * 1));
		checkRegCudaErrors(cudaMemset(d_qdiag, 0, sizeof(double) * N * 1));

		// ----------------
		// Begin regression
		// ----------------
#ifdef PRINT_CLOCKS
		printf("%d Groups, %d items per group.\n", groups, groupl);
		clock_t clock0, clock1;
		const char* time_format = "%8.2lfs";
		printf("%14s", " ");
		printf("%9s", "xtw");
		//printf("%9s", "gemm");
		printf("%9s", "Matinv");
		printf("%9s", "beta");
		printf("%9s", "S");
		printf("%9s", "C");
		printf("%9s", "betaSE");
		printf("%9s", "s_hat");
		printf("%9s", "q_diag");
		printf("%9s", "Total");
		printf("\n");
#endif
		for (size_t g = 0; g < groups; g++)
		{
#ifdef PRINT_CLOCKS
			printf("Group %6d: ", g);
			clock0 = clock();
			clock1 = clock0;
			clock_t clock_xtw = 0, clock_gemm_xtw = 0;
#endif // PRINT_CLOCKS
			size_t begin = g * groupl, end = g == (groups - 1) ? n : (g + 1) * groupl, groupn = end - begin;
			double alpha = 1.0, beta = 0.0;
			// ~~~~~~~~~~~~~~~~~~~~
			// Calculate xtwx, xtwy
			// ~~~~~~~~~~~~~~~~~~~~
			for (int i = begin; i < end; i++)
			{
				size_t e = i - begin;
#ifdef PRINT_CLOCKS
				//clock_t clocki = clock();
#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// Calculate dist, weight, xtw
				// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
				if (dm_given) { checkRegCudaErrors(cudaMemcpy(d_dists, (void *)dMat.col(i).colmem, sizeof(double) * N, cudaMemcpyHostToDevice)); }
				else { checkRegCudaErrors(gw_dist_cuda(d_dp, d_rp, N, n, i, p, theta, longlat, rp_given, d_dists, maxThreads)); }
				checkRegCudaErrors(gw_weight_cuda(bw, kernel, adaptive, d_dists, d_weights, N, 1, maxThreads));
				//checkCudaErrors(cudaDeviceSynchronize());
//#ifdef PRINT_CLOCKS
//				clock_xtw += (clock() - clocki); clocki = clock();
//#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~~~~~~~
				// Calculate xtwx, xtwy
				// ~~~~~~~~~~~~~~~~~~~~
				checkRegCudaErrors(cublasDdgmm(handle, CUBLAS_SIDE_RIGHT, K, N, d_x, K, d_weights, 1, p_xtwA[e], K));
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, K, K, N, &alpha, p_xtwA[e], K, d_x, K, &beta, p_xtwxA[e], K));
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K, 1, N, &alpha, p_xtwA[e], K, d_y, N, &beta, p_xtwyA[e], K));
#ifdef PRINT_CLOCKS
				//clock_gemm_xtw += (clock() - clocki);
#endif // PRINT_CLOCKS
			}
#ifdef PRINT_CLOCKS
			checkRegCudaErrors(cudaDeviceSynchronize());
			printf(time_format, (double)(clock() - clock1) / CLOCKS_PER_SEC); clock1 = clock();
			//printf(time_format, (double)(clock_xtw) / CLOCKS_PER_SEC);
			//printf(time_format, (double)(clock_gemm_xtw) / CLOCKS_PER_SEC);
			clock1 = clock();
#endif // PRINT_CLOCKS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~
			// Calculate inverse of xtwx
			// ~~~~~~~~~~~~~~~~~~~~~~~~~
			checkRegCudaErrors(cublasDmatinvBatched(handle, K, d_xtwxA, K, d_xtwxRA, K, d_info, groupn));
			checkRegCudaErrors(cudaDeviceSynchronize());
#ifdef PRINT_CLOCKS
			printf(time_format, (double)(clock() - clock1) / CLOCKS_PER_SEC); clock1 = clock();
#endif // PRINT_CLOCKS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Calculate Diagnosis Information
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			clock_t clock_beta = 0, clock_gemm_S = 0, clock_gemm_C = 0, clock_betasSE = 0, clock_s_hat = 0, clock_qdiag = 0;
			for (size_t i = begin; i < end; i++)
			{
				size_t e = i - begin;
#ifdef PRINT_CLOCKS
				clock_t clocki = clock();
#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~
				// Calculate beta
				// ~~~~~~~~~~~~~~
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K, 1, K, &alpha, p_xtwxRA[e], K, p_xtwyA[e], K, &beta, p_betaA[e], K));
				checkRegCudaErrors(cudaMemcpy((void *)betas.colptr(i), p_betaA[e], sizeof(double) * K, cudaMemcpyDeviceToHost));
#ifdef PRINT_CLOCKS
				clock_beta += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~~~
				// calculate ci, si
				// ~~~~~~~~~~~~~~~~
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_T, N, K, K, &alpha, p_xtwA[e], K, p_xtwxRA[e], K, &beta, d_C, N));
				checkRegCudaErrors(cudaDeviceSynchronize());
#ifdef PRINT_CLOCKS
				clock_gemm_S += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, 1, K, &alpha, d_C, N, p_xiA[i], K, &beta, d_S, N));
				checkRegCudaErrors(cudaDeviceSynchronize());
#ifdef PRINT_CLOCKS
				clock_gemm_C += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~~~~
				// calculate betasSE
				// ~~~~~~~~~~~~~~~~~
				double *betaSE = betasSE.colptr(i);
				for (size_t j = 0; j < K; j++)
				{
					checkRegCudaErrors(cublasDdot(handle, N, d_C + j * N, 1, d_C + j * N, 1, betaSE + K));
				}
#ifdef PRINT_CLOCKS
				clock_betasSE += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~~
				// calculate s_hat
				// ~~~~~~~~~~~~~~~
				double s1, s2;
				checkRegCudaErrors(cudaMemcpy((void *)&s1, d_S + i, sizeof(double), cudaMemcpyDeviceToHost));
				checkRegCudaErrors(cublasDdot(handle, N, d_S, 1, d_S, 1, &s2));
				s_hat(0) += s1;
				s_hat(1) += s2;
#ifdef PRINT_CLOCKS
				clock_s_hat += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
				// ~~~~~~~~~~~~~~~~~
				// calculate diag(q)
				// ~~~~~~~~~~~~~~~~~
				checkRegCudaErrors(gw_qdiag_cuda(d_S, N, i, d_qdiag, maxThreads));
				checkRegCudaErrors(cudaDeviceSynchronize());
#ifdef PRINT_CLOCKS
				clock_qdiag += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
			}
#ifdef PRINT_CLOCKS
			printf(time_format, (double)(clock_beta) / CLOCKS_PER_SEC);
			printf(time_format, (double)(clock_gemm_S) / CLOCKS_PER_SEC);
			printf(time_format, (double)(clock_gemm_C) / CLOCKS_PER_SEC);
			printf(time_format, (double)(clock_betasSE) / CLOCKS_PER_SEC);
			printf(time_format, (double)(clock_s_hat) / CLOCKS_PER_SEC);
			printf(time_format, (double)(clock_qdiag) / CLOCKS_PER_SEC);

			printf(time_format, (double)(clock() - clock0) / CLOCKS_PER_SEC);
			printf("\n");
#endif // PRINT_CLOCKS
		}
		checkRegCudaErrors(cudaMemcpy((double *)qdiag.mem, d_qdiag, sizeof(double) * N * 1, cudaMemcpyDeviceToHost));

		// ----------------
		// Free Cuda Memory
		// ----------------
		checkRegCudaErrors(cudaFree(d_x));
		checkRegCudaErrors(cudaFree(d_y));
		checkRegCudaErrors(cudaFree(d_dists));
		checkRegCudaErrors(cudaFree(d_weights));
		checkRegCudaErrors(cudaFree(d_xiA)); 
		checkRegCudaErrors(cudaFree(d_xtwxA)); 
		checkRegCudaErrors(cudaFree(d_xtwxRA));
		checkRegCudaErrors(cudaFree(d_info));
		checkRegCudaErrors(cudaFree(d_qdiag));
		checkRegCudaErrors(cudaFree(d_S));
		checkRegCudaErrors(cudaFree(d_C));
		delete p_xA; p_xA = nullptr;
		delete p_yA; p_yA = nullptr;
		delete p_xiA; p_xiA = nullptr;
		for (size_t i = 0; i < groupl; i++)
		{
			checkRegCudaErrors(cudaFree((void *)p_xtwA[i]));
			checkRegCudaErrors(cudaFree((void *)p_xtwxA[i]));
			checkRegCudaErrors(cudaFree((void *)p_xtwyA[i]));
			checkRegCudaErrors(cudaFree((void *)p_xtwxRA[i]));
			checkRegCudaErrors(cudaFree((void *)p_betaA[i]));
		}
		delete p_xtwA; p_xtwA = nullptr;
		delete p_xtwyA; p_xtwyA = nullptr;
		delete p_xtwxA; p_xtwxA = nullptr;
		delete p_xtwxRA; p_xtwxRA = nullptr;
		delete p_betaA; p_betaA = nullptr;
	}
	// =============================
	// Regression with out hatmatrix
	// =============================
	else
	{
		// -----------------
		// CUDA memory alloc
		// -----------------
		int groupl = smNum * maxThreads, groups = n / groupl + (n % groupl == 0 ? 0 : 1);
#ifdef PRINT_CLOCKS
		printf("%d Groups, %d items per group.\n", groups, groupl);
#endif // PRINT_CLOCKS

		// Single matrix
		double *d_x, *d_y, *d_xtw, *d_dists, *d_weights;
		checkRegCudaErrors(cudaMalloc((void **)&d_x, sizeof(double) * N * K));
		checkRegCudaErrors(cudaMalloc((void **)&d_y, sizeof(double) * N * 1));
		checkRegCudaErrors(cudaMalloc((void **)&d_xtw, sizeof(double) * N * K));
		checkRegCudaErrors(cudaMalloc((void **)&d_dists, sizeof(double) * N));
		checkRegCudaErrors(cudaMalloc((void **)&d_weights, sizeof(double) * N));
		checkRegCudaErrors(cudaMemcpy(d_x, x.mem, sizeof(double) * N * K, cudaMemcpyHostToDevice));
		checkRegCudaErrors(cudaMemcpy(d_y, y.mem, sizeof(double) * N * 1, cudaMemcpyHostToDevice));
		double **p_xtwxA, **p_xtwyA, **p_xtwxRA, **p_betaA;
		p_xtwxA = new double*[groupl];
		p_xtwyA = new double*[groupl];
		p_xtwxRA = new double*[groupl];
		p_betaA = new double*[groupl];
		for (size_t i = 0; i < groupl; i++)
		{
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwxA[i], sizeof(double) * K * K));
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwyA[i], sizeof(double) * K * 1));
			checkRegCudaErrors(cudaMalloc((void **)&p_xtwxRA[i], sizeof(double) * K * K));
			checkRegCudaErrors(cudaMalloc((void **)&p_betaA[i], sizeof(double) * K * 1));
		}
		double **d_xtwxA, **d_xtwxRA;
		checkRegCudaErrors(cudaMalloc((void **)&d_xtwxA, sizeof(double*) * groupl));
		checkRegCudaErrors(cudaMalloc((void **)&d_xtwxRA, sizeof(double*) * groupl));
		checkRegCudaErrors(cudaMemcpy(d_xtwxA, p_xtwxA, sizeof(double*) * groupl, cudaMemcpyHostToDevice));
		checkRegCudaErrors(cudaMemcpy(d_xtwxRA, p_xtwxRA, sizeof(double*) * groupl, cudaMemcpyHostToDevice));
		int *d_info;
		checkRegCudaErrors(cudaMalloc((void **)&d_info, sizeof(int) * groupl));

#ifdef PRINT_CLOCKS
		// clocks
		// ------
		clock_t clock0 = clock();
		const char* time_format = "%10.2lfs";
		printf("%13s", " ");
		printf("%11s", "xtw");
		printf("%11s", "gemm");
		printf("%11s", "Matinv");
		printf("%11s", "beta");
		printf("%11s", "Total");
		printf("\n");
#endif // PRINT_CLOCKS

		// ----------------
		// Begin regression
		// ----------------
		//printf("\n");
		double one = 1.0, zero = 0.0;
		for (size_t g = 0; g < groups; g++)
		{
			size_t begin = g * groupl, end = g == (groups - 1) ? n : (g + 1) * groupl, groupn = end - begin;
			// ~~~~~~~~~~~~~~~~~~~
			// Calculate xtwx xtwy
			// ~~~~~~~~~~~~~~~~~~~
#ifdef PRINT_CLOCKS
			printf("Group %5d: ", g);
			clock_t clock1 = clock(), clock_gemm = 0;
#endif
			for (int i = begin; i < end; i++)
			{
				size_t e = i - begin;
#ifdef PRINT_CLOCKS
				clock_t clocki = clock();
#endif // PRINT_CLOCKS
				if (dm_given) { checkRegCudaErrors(cudaMemcpy(d_dists, (void *)dMat.col(i).colmem, sizeof(double) * N, cudaMemcpyHostToDevice)); }
				else { checkRegCudaErrors(gw_dist_cuda(d_dp, d_rp, N, n, i, p, theta, longlat, rp_given, d_dists, maxThreads)); }
				checkRegCudaErrors(gw_weight_cuda(bw, kernel, adaptive, d_dists, d_weights, N, 1, maxThreads));
				checkRegCudaErrors(gw_xtw_cuda(d_x, d_weights, N, K, d_xtw, maxThreads));
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, K, K, N, &one, d_xtw, K, d_x, K, &zero, p_xtwxA[e], K));
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K, 1, N, &one, d_xtw, K, d_y, N, &zero, p_xtwyA[e], K));
#ifdef PRINT_CLOCKS
				clock_gemm += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
			}
#ifdef PRINT_CLOCKS
			printf(time_format, (double)(clock_gemm) / CLOCKS_PER_SEC); clock1 = clock();
#endif
			// ~~~~~~~~~~~~~~
			// Matrix inverse
			// ~~~~~~~~~~~~~~
			checkRegCudaErrors(cublasDmatinvBatched(handle, K, d_xtwxA, K, d_xtwxRA, K, d_info, groupn));
#ifdef PRINT_CLOCKS
			printf(time_format, (double)(clock() - clock1) / CLOCKS_PER_SEC); clock1 = clock();
#endif
			// ~~~~~~~~~~~~~~
			// calculate beta
			// ~~~~~~~~~~~~~~
			for (size_t i = begin; i < end; i++)
			{
				size_t e = i - begin;
				checkRegCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K, 1, K, &one, p_xtwxRA[e], K, p_xtwyA[e], K, &zero, p_betaA[e], K));
				checkRegCudaErrors(cudaMemcpy((void *)betas.col(i).colmem, p_betaA[e], sizeof(double) * K, cudaMemcpyDeviceToHost));
			}
#ifdef PRINT_CLOCKS
			printf(time_format, (double)(clock() - clock1) / CLOCKS_PER_SEC); clock1 = clock();
			printf(time_format, (double)(clock() - clock0) / CLOCKS_PER_SEC); clock0 = clock();
			printf("\n");
#endif
		}
#ifdef PRINT_CLOCKS
		printf(time_format, (double)(clock() - clock0) / CLOCKS_PER_SEC);
		printf("\n");
#endif
		// ----------------
		// Free Cuda Memory
		// ----------------
		checkRegCudaErrors(cudaFree(d_x));
		checkRegCudaErrors(cudaFree(d_y));
		checkRegCudaErrors(cudaFree(d_dists));
		checkRegCudaErrors(cudaFree(d_weights));
		checkRegCudaErrors(cudaFree(d_xtwxA));
		checkRegCudaErrors(cudaFree(d_xtwxRA));
		checkRegCudaErrors(cudaFree(d_info));
		for (size_t i = 0; i < groupl; i++)
		{
			checkRegCudaErrors(cudaFree(p_xtwxA[i]));
			checkRegCudaErrors(cudaFree(p_xtwyA[i]));
			checkRegCudaErrors(cudaFree(p_xtwxRA[i]));
			checkRegCudaErrors(cudaFree(p_betaA[i]));
		}
		delete[] p_xtwxA; p_xtwxA = nullptr;
		delete[] p_xtwyA; p_xtwyA = nullptr;
		delete[] p_xtwxRA; p_xtwxRA = nullptr;
		delete[] p_betaA; p_betaA = nullptr;
	}
	checkRegCudaErrors(cudaFree(d_dp));
	if (rp_given)
	{
		checkRegCudaErrors(cudaFree(d_rp));
	}
	cublasDestroy(handle);
	return true;
}

double CGWmodelCUDA::CV(double p, double theta, bool longlat, double bw, int kernel, bool adaptive, int groupl, int gpuID)
{
	int n = rp.n_rows, N = dp.n_rows, K = x.n_rows;
	// ====================
	// Check GPU properties
	// ====================
	int gpuCount = 0;
	checkCvCudaErrors(cudaGetDeviceCount(&gpuCount));
	if (gpuID >= gpuCount)
	{
		return false;
	}
	cudaDeviceProp devProp;
	checkCvCudaErrors(cudaGetDeviceProperties(&devProp, gpuID));
	int smNum = devProp.multiProcessorCount;
	int maxThreads = devProp.maxThreadsPerBlock;
	if (groupl <= 0)
	{
		groupl = smNum * maxThreads / 2;
	}

#ifdef PRINT_CLOCKS
	printf("GPU Device: %s. \t Group length: %d\n", devProp.name, groupl);
#endif // PRINT_CLOCKS


	// ============
	// Cuda prepare
	// ============
	cublasHandle_t handle;
	cublasCreate(&handle);

	// ====================
	// Armadillo Precompute
	// ====================
	//betas = mat(K, n, fill::zeros);
	double cv = 0.0;

	// ============
	// dp rp prepare
	// ============
	double *d_dp;
	checkCvCudaErrors(cudaMalloc((void **)&d_dp, sizeof(double) * N * 2));
	checkCvCudaErrors(cudaMemcpy(d_dp, dp.mem, sizeof(double) * N * 2, cudaMemcpyHostToDevice));
	if (p != 2.0 && theta != 0.0 && !longlat)
	{
		checkCvCudaErrors(gw_coordinate_rotate_cuda(d_dp, N, theta, maxThreads));
	}
	int groups = n / groupl + (n % groupl == 0 ? 0 : 1);
#ifdef PRINT_CLOCKS
	printf("%d Groups, %d items per group.\n", groups, groupl);
#endif // PRINT_CLOCKS

	// Single matrix
	double *d_x, *d_y, *d_xtw, *d_dists, *d_weights;
	checkCvCudaErrors(cudaMalloc((void **)&d_x, sizeof(double) * N * K));
	checkCvCudaErrors(cudaMalloc((void **)&d_y, sizeof(double) * N * 1));
	checkCvCudaErrors(cudaMalloc((void **)&d_xtw, sizeof(double) * N * K));
	checkCvCudaErrors(cudaMalloc((void **)&d_dists, sizeof(double) * N));
	checkCvCudaErrors(cudaMalloc((void **)&d_weights, sizeof(double) * N));
	checkCvCudaErrors(cudaMemcpy(d_x, x.mem, sizeof(double) * N * K, cudaMemcpyHostToDevice));
	checkCvCudaErrors(cudaMemcpy(d_y, y.mem, sizeof(double) * N * 1, cudaMemcpyHostToDevice));
	double **p_xtwxA, **p_xtwyA, **p_xtwxRA, **p_betaA;
	p_xtwxA = new double*[groupl];
	p_xtwyA = new double*[groupl];
	p_xtwxRA = new double*[groupl];
	p_betaA = new double*[groupl];
	for (size_t i = 0; i < groupl; i++)
	{
		checkCvCudaErrors(cudaMalloc((void **)&p_xtwxA[i], sizeof(double) * K * K));
		checkCvCudaErrors(cudaMalloc((void **)&p_xtwyA[i], sizeof(double) * K * 1));
		checkCvCudaErrors(cudaMalloc((void **)&p_xtwxRA[i], sizeof(double) * K * K));
		checkCvCudaErrors(cudaMalloc((void **)&p_betaA[i], sizeof(double) * K * 1));
	}
	double **p_xiA = new double*[N];
	for (size_t i = 0; i < N; i++)
	{
		p_xiA[i] = d_x + i * K;
	}
	double **d_xtwxA, **d_xtwxRA;
	checkCvCudaErrors(cudaMalloc((void **)&d_xtwxA, sizeof(double*) * groupl));
	checkCvCudaErrors(cudaMalloc((void **)&d_xtwxRA, sizeof(double*) * groupl));
	checkCvCudaErrors(cudaMemcpy(d_xtwxA, p_xtwxA, sizeof(double*) * groupl, cudaMemcpyHostToDevice));
	checkCvCudaErrors(cudaMemcpy(d_xtwxRA, p_xtwxRA, sizeof(double*) * groupl, cudaMemcpyHostToDevice));
	int *p_info, *d_info;
	p_info = new int[groupl];
	checkCvCudaErrors(cudaMalloc((void **)&d_info, sizeof(int) * groupl));

#ifdef PRINT_CLOCKS
	// clocks
	// ------
	clock_t clock0 = clock();
	const char* time_format = "%10.2lfs";
	printf("%13s", " ");
	printf("%11s", "xtw");
	printf("%11s", "gemm");
	printf("%11s", "Matinv");
	printf("%11s", "beta");
	printf("%11s", "Total");
	printf("\n");
#endif // PRINT_CLOCKS

	// ----------------
	// Begin regression
	// ----------------
	//printf("\n");
	double one = 1.0, zero = 0.0;
	for (size_t g = 0; g < groups; g++)
	{
		size_t begin = g * groupl, end = g == (groups - 1) ? n : (g + 1) * groupl, groupn = end - begin;
		// ~~~~~~~~~~~~~~~~~~~
		// Calculate xtwx xtwy
		// ~~~~~~~~~~~~~~~~~~~
#ifdef PRINT_CLOCKS
		printf("Group %5d: ", g);
		clock_t clock1 = clock(), clock_gemm = 0;
#endif
		for (int i = begin; i < end; i++)
		{
			size_t e = i - begin;
#ifdef PRINT_CLOCKS
			clock_t clocki = clock();
#endif // PRINT_CLOCKS
			if (dm_given) { checkCvCudaErrors(cudaMemcpy(d_dists, (void *)dMat.col(i).colmem, sizeof(double) * N, cudaMemcpyHostToDevice)); }
			else { checkCvCudaErrors(gw_dist_cuda(d_dp, d_dp, N, n, i, p, theta, longlat, rp_given, d_dists, maxThreads)); }
			checkCvCudaErrors(gw_weight_cuda(bw, kernel, adaptive, d_dists, d_weights, N, 1, maxThreads));
			checkCvCudaErrors(cudaMemcpy(d_weights + i, &zero, sizeof(double), cudaMemcpyHostToDevice));
			checkCvCudaErrors(gw_xtw_cuda(d_x, d_weights, N, K, d_xtw, maxThreads));
			checkCvCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, K, K, N, &one, d_xtw, K, d_x, K, &zero, p_xtwxA[e], K));
			checkCvCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K, 1, N, &one, d_xtw, K, d_y, N, &zero, p_xtwyA[e], K));
#ifdef PRINT_CLOCKS
			clock_gemm += (clock() - clocki); clocki = clock();
#endif // PRINT_CLOCKS
		}
#ifdef PRINT_CLOCKS
		printf(time_format, (double)(clock_gemm) / CLOCKS_PER_SEC); clock1 = clock();
#endif
		// ~~~~~~~~~~~~~~
		// Matrix inverse
		// ~~~~~~~~~~~~~~
		checkCvCudaErrors(cublasDmatinvBatched(handle, K, d_xtwxA, K, d_xtwxRA, K, d_info, groupn));
		checkCvCudaErrors(cudaMemcpy(p_info, d_info, sizeof(int) * groupl, cudaMemcpyDeviceToHost));
		for (size_t i = begin; i < end; i++)
		{
			size_t e = i - begin;
			int info = p_info[e];
			if (info > 0)
			{
				#ifndef WIN32
				throw std::exception(std::runtime_error("Matrix seems singular."));
				#else
				throw std::exception("Matrix seems singular.");
				#endif
			}
		}
#ifdef PRINT_CLOCKS
		printf(time_format, (double)(clock() - clock1) / CLOCKS_PER_SEC); clock1 = clock();
#endif
		// ~~~~~~~~~~~~~~
		// calculate beta
		// ~~~~~~~~~~~~~~
		for (size_t i = begin; i < end; i++)
		{
			size_t e = i - begin;
			double yhat = 0;
			checkCvCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K, 1, K, &one, p_xtwxRA[e], K, p_xtwyA[e], K, &zero, p_betaA[e], K));
			checkCvCudaErrors(cublasDdot(handle, K, p_xiA[i], 1, p_betaA[e], 1, &yhat));
			cv += ((y(i) - yhat) * (y(i) - yhat));
		}
#ifdef PRINT_CLOCKS
		printf(time_format, (double)(clock() - clock1) / CLOCKS_PER_SEC); clock1 = clock();
		printf(time_format, (double)(clock() - clock0) / CLOCKS_PER_SEC); clock0 = clock();
		printf("\n");
#endif
	}
#ifdef PRINT_CLOCKS
	printf(time_format, (double)(clock() - clock0) / CLOCKS_PER_SEC);
	printf("\n");
#endif
	// ----------------
	// Free Cuda Memory
	// ----------------
	checkCvCudaErrors(cudaFree(d_x));
	checkCvCudaErrors(cudaFree(d_y));
	checkCvCudaErrors(cudaFree(d_dists));
	checkCvCudaErrors(cudaFree(d_weights));
	checkCvCudaErrors(cudaFree(d_xtwxA));
	checkCvCudaErrors(cudaFree(d_xtwxRA));
	checkCvCudaErrors(cudaFree(d_info));
	for (size_t i = 0; i < groupl; i++)
	{
		checkCvCudaErrors(cudaFree(p_xtwxA[i]));
		checkCvCudaErrors(cudaFree(p_xtwyA[i]));
		checkCvCudaErrors(cudaFree(p_xtwxRA[i]));
		checkCvCudaErrors(cudaFree(p_betaA[i]));
	}
	delete[] p_xtwxA; p_xtwxA = nullptr;
	delete[] p_xtwxRA; p_xtwxRA = nullptr;
	delete[] p_xtwyA; p_xtwyA = nullptr;
	delete[] p_betaA; p_betaA = nullptr;
	checkCvCudaErrors(cudaFree(d_dp));
	cublasDestroy(handle);
	return cv;
}
