#define M_PI       3.14159265358979323846
#define DOUBLE_EPS 1e-12

#include "GWmodelKernel.h"
#include <device_launch_parameters.h>
#include <cmath>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#define POWDI(x,i) pow(x,i)
#define GAUSSIAN 0
#define EXPONENTIAL 1
#define BISQUARE 2
#define TRICUBE 3
#define BOXCAR 4

__global__ void coordinate_rotate(double* coords, int n, double costheta, double sintheta)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double ix = coords[i], iy = coords[i + n];
	double ox = ix * costheta - iy * sintheta;
	double oy = ix * sintheta + iy * costheta;
	coords[i] = ox;
	coords[i + n] = oy;
}

cudaError_t gw_coordinate_rotate_cuda(double* d_coords, int n, double theta, int threads)
{
	cudaError_t error;
	dim3 blockSize(threads), gridSize((n + blockSize.x - 1) / blockSize.x);
	coordinate_rotate << <gridSize, blockSize >> > (d_coords, n, cos(theta), sin(theta));
	error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		return error;
	}
	return cudaSuccess;
}

__global__ void eu_dist_vec_kernel(const double *dp, int ndp, const double *rp, int focus, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int i = index;
	double ix = *(dp + i), iy = *(dp + i + ndp);
	double ox = *(rp + focus), oy = *(rp + focus + nrp);
	double dist = hypot((ix - ox), (iy - oy));
	*(dists + i) = dist;
}

__global__ void cd_dist_vec_kernel(const double *dp, int ndp, const double *rp, int focus, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int i = index;
	double ix = dp[i], iy = dp[i + ndp];
	double ox = *(rp + focus), oy = *(rp + focus + nrp);
	double dist = fmax(fabs(ix - ox), fabs(iy - oy));
	*(dists + i) = dist;
}

__global__ void md_dist_vec_kernel(const double *dp, int ndp, const double *rp, int focus, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int i = index;
	double ix = dp[i], iy = dp[i + ndp];
	double ox = *(rp + focus), oy = *(rp + focus + nrp);
	double dist = fabs(ix - ox) + fabs(iy - oy);
	*(dists + i) = dist;
}

__global__ void mk_dist_vec_kernel(const double *dp, int ndp, const double *rp, int focus, int nrp, double p, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int i = index;
	double ix = dp[i], iy = dp[i + ndp];
	double ox = *(rp + focus), oy = *(rp + focus + nrp);
	double dist = pow(pow(fabs(ix - ox), p) + pow(fabs(iy - oy), p), 1.0 / p);
	*(dists + i) = dist;
}

__global__ void eu_dist_mat_kernel(const double *dp, int ndp, const double *rp, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int r = index / nrp, d = index % ndp;
	double ix = dp[d], iy = dp[d + ndp], ox = rp[r], oy = rp[r + nrp];
	double dist = hypot((ix - ox), (iy - oy));
	*(dists + index) = dist;
}

__global__ void cd_dist_mat_kernel(const double *dp, int ndp, const double *rp, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int r = index / nrp, d = index % ndp;
	double ix = dp[d], iy = dp[d + ndp], ox = rp[r], oy = rp[r + nrp];
	double dist = fmax(fabs(ix - ox), fabs(iy - oy));
	*(dists + index) = dist;
}

__global__ void md_dist_mat_kernel(const double *dp, int ndp, const double *rp, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int r = index / nrp, d = index % ndp;
	double ix = dp[d], iy = dp[d + ndp], ox = rp[r], oy = rp[r + nrp];
	double dist = fabs(ix - ox) + fabs(iy - oy);
	*(dists + index) = dist;
}

__global__ void mk_dist_mat_kernel(const double *dp, int ndp, const double *rp, int nrp, double p, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int r = index / nrp, d = index % ndp;
	double ix = dp[d], iy = dp[d + ndp], ox = rp[r], oy = rp[r + nrp];
	double dist = pow(pow(fabs(ix - ox), p) + pow(fabs(iy - oy), p), 1.0 / p);
	*(dists + index) = dist;
}

__device__ double sp_gcdist(double lon1, double lon2, double lat1, double lat2)
{

	double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
	double w, R, a, f, D, H1, H2;
	double lat1R, lat2R, lon1R, lon2R, DE2RA;

	DE2RA = M_PI / 180;
	a = 6378.137;              /* WGS-84 equatorial radius in km */
	f = 1.0 / 298.257223563;     /* WGS-84 ellipsoid flattening factor */

	if (fabs(lat1 - lat2) < DOUBLE_EPS)
	{
		if (fabs(lon1 - lon2) < DOUBLE_EPS)
		{
			return 0.0;
			/* Wouter Buytaert bug caught 100211 */
		}
		else if (fabs((fabs(lon1) + fabs(lon2)) - 360.0) < DOUBLE_EPS)
		{
			return 0.0;
		}
	}
	lat1R = lat1 * DE2RA;
	lat2R = lat2 * DE2RA;
	lon1R = lon1 * DE2RA;
	lon2R = lon2 * DE2RA;

	F = (lat1R + lat2R) / 2.0;
	G = (lat1R - lat2R) / 2.0;
	L = (lon1R - lon2R) / 2.0;

	sinG2 = POWDI(sin(G), 2);
	cosG2 = POWDI(cos(G), 2);
	sinF2 = POWDI(sin(F), 2);
	cosF2 = POWDI(cos(F), 2);
	sinL2 = POWDI(sin(L), 2);
	cosL2 = POWDI(cos(L), 2);

	S = sinG2 * cosL2 + cosF2 * sinL2;
	C = cosG2 * cosL2 + sinF2 * sinL2;

	w = atan(sqrt(S / C));
	R = sqrt(S*C) / w;

	D = 2 * w*a;
	H1 = (3 * R - 1) / (2 * C);
	H2 = (3 * R + 1) / (2 * S);

	return D * (1 + f * H1*sinF2*cosG2 - f * H2*cosF2*sinG2);
}

__global__ void sp_dist_vec_kernel(const double *dp, int ndp, const double *rp, int focus, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int i = index;
	double ix = dp[i], iy = dp[i + ndp];
	double ox = *(rp + focus), oy = *(rp + focus + nrp);
	dists[i] = sp_gcdist(ix, ox, iy, oy);
}

__global__ void sp_dist_mat_kernel(const double *dp, int ndp, const double *rp, int nrp, double *dists)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= ndp) return;
	int r = index / nrp, d = index % ndp;
	double ix = dp[d], iy = dp[d + ndp], ox = rp[r], oy = rp[r + nrp];
	dists[index] = sp_gcdist(ix, ox, iy, oy);
}

cudaError_t gw_dist_cuda(double *d_dp, double *d_rp, int ndp, int nrp, int focus, double p, double theta, bool longlat, bool rp_given, double *d_dists, int threads)
{
	cudaError_t error;
	int isFocus = focus > -1;
	if (isFocus)
	{
		dim3 blockSize(threads), gridSize((ndp + blockSize.x - 1) / blockSize.x);
		if (longlat)
		{
			sp_dist_vec_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, focus, nrp, d_dists);
		}
		else
		{
			if (p == 2.0)
				eu_dist_vec_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, focus, nrp, d_dists);
			else if (p == 1.0)
				cd_dist_vec_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, focus, nrp, d_dists);
			else if (p == -1.0)
				md_dist_vec_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, focus, nrp, d_dists);
			else
				mk_dist_vec_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, focus, nrp, p, d_dists);
		}
	}
	else
	{
		if (rp_given)
		{
			dim3 blockSize(threads), gridSize((ndp * nrp + blockSize.x - 1) / blockSize.x);
			if (longlat)
				sp_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, nrp, d_dists);
			else
			{
				if (p == 2.0)
					eu_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, nrp, d_dists);
				else if (p == 1.0)
					cd_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, nrp, d_dists);
				else if (p == -1.0)
					md_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, nrp, d_dists);
				else
					mk_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_rp, nrp, p, d_dists);
			}
		}
		else
		{
			dim3 blockSize(threads), gridSize((ndp * ndp + blockSize.x - 1) / blockSize.x);
			if (longlat)
				sp_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_dp, ndp, d_dists);
			else
			{
				if (p == 2.0)
					eu_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_dp, ndp, d_dists);
				else if (p == 1.0)
					cd_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_dp, ndp, d_dists);
				else if (p == -1.0)
					md_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_dp, ndp, d_dists);
				else
					mk_dist_mat_kernel << <gridSize, blockSize >> > (d_dp, ndp, d_dp, ndp, p, d_dists);
			}
		}
	}
	error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		return error;
	}
	return cudaSuccess;
}

__global__ void gw_weight_gaussian_kernel(const double *d_dists, double bw, double *d_weights, int n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double dist = d_dists[i];
	d_weights[i] = exp((dist * dist) / ((-2)*(bw * bw)));
}

__global__ void gw_weight_exponential_kernel(const double *d_dists, double bw, double *d_weights, int n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double dist = d_dists[i];
	d_weights[i] = exp(-dist / bw);
}

__global__ void gw_weight_bisquare_kernel(const double *d_dists, double bw, double *d_weights, int n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double dist = d_dists[i];
	d_weights[i] = dist > bw ? 0 : (1 - (dist * dist) / (bw * bw))*(1 - (dist * dist) / (bw * bw));
}

__global__ void gw_weight_tricube_kernel(const double *d_dists, double bw, double *d_weights, int n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double dist = d_dists[i];
	d_weights[i] = dist > bw ? 0 :
		(1 - (dist * dist * dist) / (bw * bw * bw))*
		(1 - (dist * dist * dist) / (bw * bw * bw))*
		(1 - (dist * dist * dist) / (bw * bw * bw));
}

__global__ void gw_weight_boxcar_kernel(const double *d_dists, double bw, double *d_weights, int n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double dist = d_dists[i];
	d_weights[i] = dist > bw ? 0 : 1;
}

typedef void(*WEIGHT_KERNEL_CUDA)(const double*, double, double*, int);

const WEIGHT_KERNEL_CUDA GWRKernelCuda[5] = {
  gw_weight_gaussian_kernel,
  gw_weight_exponential_kernel,
  gw_weight_bisquare_kernel,
  gw_weight_tricube_kernel,
  gw_weight_boxcar_kernel
};

cudaError_t gw_weight_cuda(double bw, int kernel, bool adaptive, double *d_dists, double *d_weight, int ndp, int nrp, int threads)
{
	cudaError_t error;
	const WEIGHT_KERNEL_CUDA *kerf = GWRKernelCuda + kernel;
	switch (adaptive)
	{
		case true:
		{
			dim3 blockSize(threads), gridSize((ndp + blockSize.x - 1) / blockSize.x);
			for (size_t f = 0; f < nrp; f++)
			{
				// Backup d_dists, used for sort
				double *d_dists_bak;
				cudaMalloc((void **)&d_dists_bak, sizeof(double) * ndp);
				cudaMemcpy(d_dists_bak, d_dists + f * ndp, sizeof(double) * ndp, cudaMemcpyDeviceToDevice);
				thrust::device_ptr<double> v_dists(d_dists_bak);
				thrust::sort(v_dists, v_dists + ndp);
				// Calculate weight for each distance
				double bw_dist = v_dists[(int)(bw < ndp ? bw : ndp) - 1];
				(*kerf) << <gridSize, blockSize >> > (d_dists + f * ndp, bw_dist, d_weight + f * ndp, ndp);
				// Free d_dists_bak
				cudaFree(d_dists_bak);
				d_dists_bak = nullptr;
				// Get error
				error = cudaGetLastError();
				if (error != cudaSuccess)
				{
					return error;
				}
			}
			break;
		}
		default:
		{
			dim3 blockSize(threads), gridSize((ndp * nrp + blockSize.x - 1) / blockSize.x);
			(*kerf) << <gridSize, blockSize >> > (d_dists, bw, d_weight, ndp * nrp);
			error = cudaGetLastError();
			if (error != cudaSuccess)
			{
				return error;
			}
			break;
		}
	}
	return cudaSuccess;
}

__global__ void gw_xtw_kernel(const double* d_x, const double* d_wights, int n, int k, double* d_xtw)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	double weight = d_wights[index];
	for (int j = 0; j < k; j++)
	{
		int p = j + i * k;
		d_xtw[p] = d_x[p] * weight;
	}
}

cudaError_t gw_xtw_cuda(const double* d_x, const double* d_weight, int n, int k, double* d_xtw, int threads)
{
	cudaError_t error;
	dim3 blockSize(threads), gridSize((n + blockSize.x - 1) / blockSize.x);
	gw_xtw_kernel << <gridSize, blockSize >> > (d_x, d_weight, n, k, d_xtw);
	error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		return error;
	}
	return cudaSuccess;
}

__global__ void gw_xdy_kernel(const double* x, const double* y, int n, double* xdy)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	int i = index;
	xdy[i] = x[i] * y[i];
}

cudaError_t gw_xdy_cuda(const double* d_x, const double* d_y, int n, double * d_xdoty, int threads)
{
	cudaError_t error;
	dim3 blockSize(threads), gridSize((n + blockSize.x - 1) / blockSize.x);
	gw_xdy_kernel << <gridSize, blockSize >> > (d_x, d_y, n, d_xdoty);
	error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		return error;
	}
	return cudaSuccess;
}

__global__ void gw_xdx_kernel(const double* x, int n, double* xdx)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= n) return;
	double a = x[index];
	xdx[index] = a * a;
}

cudaError_t gw_xdx_cuda(const double* d_x, int n, double * d_xdotx, int threads)
{
	cudaError_t error;
	dim3 blockSize(threads), gridSize((n + blockSize.x - 1) / blockSize.x);
	gw_xdx_kernel << <gridSize, blockSize >> > (d_x, n, d_xdotx);
	error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		return error;
	}
	return cudaSuccess;
}

__global__ void gw_qdiag_kernel(const double* d_si, int n, int p, double* d_q)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int i = ((index >= n) ? n - 1 : index);
	d_q[i] += ((i == p) ? (1 - d_si[i])*(1 - d_si[i]) : d_si[i] * d_si[i]);
}

cudaError_t gw_qdiag_cuda(const double* d_si, int n, int p, double* d_q, int threads)
{
	cudaError_t error;
	dim3 blockSize(threads), gridSize((n + blockSize.x - 1) / blockSize.x);
	gw_qdiag_kernel << <gridSize, blockSize >> > (d_si, n, p, d_q);
	error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		return error;
	}
	return cudaSuccess;
}