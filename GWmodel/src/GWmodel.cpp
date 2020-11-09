// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

#define POWDI(x,i) pow(x,i)

#define GAUSSIAN 0
#define EXPONENTIAL 1
#define BISQUARE 2
#define TRICUBE 3
#define BOXCAR 4

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------

const mat mSum(2, 1, fill::ones);

//distance matrix calculation
//coords must be a matrix with 2 columns
// [[Rcpp::export]]
mat coordinate_rotate(mat coords, double theta)
{
	int n = coords.n_rows;
	mat rotated_coords(n, 2);
	rotated_coords.col(0) = coords.col(0) * cos(theta) - coords.col(1) * sin(theta);
	rotated_coords.col(1) = coords.col(0) * sin(theta) + coords.col(1) * cos(theta);
	return rotated_coords;
}

//Eudclidean distance matrix
// [[Rcpp::export]]
mat eu_dist_mat(mat in_locs, mat out_locs)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat eu_dist(n_in, n_out);
	int i = 0, j = 0;
	for (i = 0; i < n_in; i++)
	{
		for (j = 0; j < n_out; j++)
		{
		  eu_dist(i,j) = sum(pow(in_locs.row(i) - out_locs.row(j),2));
		}
	}
	return sqrt(eu_dist);
}
//symmetrical distance matrix
// [[Rcpp::export]]
mat eu_dist_smat(mat in_locs)
{
	int n = in_locs.n_rows;
	mat eu_dist(n, n);
	for (int k = 0; k < n * n; k++)
	{
		int i = k / n, j = k % n;
		eu_dist(i, j) = sum(pow(in_locs.row(i) - in_locs.row(j), 2));
		eu_dist(j, i) = eu_dist(i, j);
	}
	return sqrt(eu_dist);
}

// [[Rcpp::export]]
vec eu_dist_vec(mat in_locs, vec out_loc)
{
	int n_in = in_locs.n_rows;
	vec eu_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		eu_dist(i) = sum(pow(in_locs.row(i) - trans(out_loc), 2));
	}
	return sqrt(eu_dist);
	// mat v_span(n_in, 1, fill::ones);
	// mat m_diff = in_locs - v_span * trans(out_loc);
	// return sqrt(m_diff % m_diff * mSum);
}

//Manhattan distance matrix
// [[Rcpp::export]]
mat md_dist_mat(mat in_locs, mat out_locs)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat md_dist(n_in, n_out);
	for (int i = 0; i < n_in; i++)
	{
		for (int j = 0; j < n_out; j++)
		{
			md_dist(i, j) = sum(abs(in_locs.row(i) - out_locs.row(j)));
		}
	}
	return md_dist;
}

//symmetrical distance matrix
// [[Rcpp::export]]
mat md_dist_smat(mat in_locs)
{
	int n = in_locs.n_rows;
	mat md_dist(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			md_dist(i, j) = sum(abs(in_locs.row(i) - in_locs.row(j)));
			md_dist(j, i) = md_dist(i, j);
		}
	}
	return md_dist;
}
// [[Rcpp::export]]
vec md_dist_vec(mat in_locs, vec out_loc)
{
	int n_in = in_locs.n_rows;
	vec md_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		md_dist(i) = sum(abs(in_locs.row(i) - trans(out_loc)));
	}
	return md_dist;
}

//Chebyshev distance matrix
// [[Rcpp::export]]
mat cd_dist_mat(mat in_locs, mat out_locs)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat cd_dist(n_in, n_out);
	for (int i = 0; i < n_in; i++)
	{
		for (int j = i; j < n_out; j++)
		{
			cd_dist(i, j) = max(abs(in_locs.row(i) - out_locs.row(j)));
			cd_dist(j, i) = cd_dist(i, j);
		}
	}
	return cd_dist;
}

//symmetrical distance matrix
// [[Rcpp::export]]
mat cd_dist_smat(mat in_locs)
{
	int n = in_locs.n_rows;
	mat cd_dist(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			cd_dist(i, j) = max(abs(in_locs.row(i) - in_locs.row(j)));
			cd_dist(j, i) = cd_dist(i, j);
		}
	}
	return cd_dist;
}
// [[Rcpp::export]]
vec cd_dist_vec(mat in_locs, vec out_loc)
{
	int n_in = in_locs.n_rows;
	vec cd_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		cd_dist(i) = max(abs(in_locs.row(i) - trans(out_loc)));
	}
	return cd_dist;
}

//Minkowski distance matrix
// [[Rcpp::export]]
mat mk_dist_mat(mat in_locs, mat out_locs, double p)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat mk_dist(n_in, n_out);
	for (int i = 0; i < n_in; i++)
	{
		for (int j = 0; j < n_out; j++)
		{
			mk_dist(i, j) = pow(sum(pow(abs(in_locs.row(i) - out_locs.row(j)), p)), 1.0 / p);
		}
	}
	return mk_dist;
}
//sqrt(sum(pow(in_locs.row(i) - trans(out_loc),2)))
//symmetrical distance matrix
// [[Rcpp::export]]
mat mk_dist_smat(mat in_locs, double p)
{
	int n = in_locs.n_rows;
	mat mk_dist(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			mk_dist(i, j) = pow(sum(pow(abs(in_locs.row(i) - in_locs.row(j)), p)), 1.0 / p);
			mk_dist(j, i) = mk_dist(i, j);
		}
	}
	return mk_dist;
}

// [[Rcpp::export]]
vec mk_dist_vec(mat in_locs, vec out_loc, double p)
{
	int n_in = in_locs.n_rows;
	vec mk_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		mk_dist(i) = pow(sum(pow(abs(in_locs.row(i) - trans(out_loc)), p)), 1.0 / p);
	}
	return mk_dist;
}
//Weight matrix
//Bisuqare weight
// [[Rcpp::export]]
vec bisq_wt_vec(vec distv, double bw)
{
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for (int i = 0; i < n; i++)
	{
		if (distv(i) <= bw)
			wtv(i) = pow(1 - pow(distv(i) / bw, 2), 2);
	}
	return wtv;
}
//Calculated by column, the length of bw must be the same the number of columns of distm
// [[Rcpp::export]]
mat bisq_wt_mat(mat distm, vec bw)
{
	int m = distm.n_cols;
	int n = distm.n_rows;
	mat wtm;
	wtm.zeros(n, m);

	for (int k = 0; k < m * n; k++)
	{
		int i = k / n, j = k % n;
		if (distm(j, i) <= bw(i))
			wtm(j, i) = (1 - distm(j, i) / bw(i) * distm(j, i) / bw(i)) * (1 - distm(j, i) / bw(i) * distm(j, i) / bw(i));
	}
	return wtm;
}

//Gaussian weight
// [[Rcpp::export]]
vec gauss_wt_vec(vec distv, double bw)
{
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for (int i = 0; i < n; i++)
	{
		wtv(i) = exp(pow(distv(i), 2) / ((-2) * pow(bw, 2)));
	}
	return wtv;
}
// [[Rcpp::export]]
mat gauss_wt_mat(mat distm, vec bw)
{
	int m = distm.n_cols;
	int n = distm.n_rows;
	mat wtm;
	wtm.zeros(n, m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			wtm(j, i) = exp(pow(distm(j, i), 2) / ((-2) * pow(bw(i), 2)));
		}
	}
	return wtm;
}

//Tricube weight
// [[Rcpp::export]]
vec tri_wt_vec(vec distv, double bw)
{
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for (int i = 0; i < n; i++)
	{
		if (distv(i) <= bw)
			wtv(i) = pow(1 - pow(distv(i), 3) / pow(bw, 3), 3);
	}
	return wtv;
}
// [[Rcpp::export]]
mat tri_wt_mat(mat distm, vec bw)
{
	int m = distm.n_cols;
	int n = distm.n_rows;
	mat wtm;
	wtm.zeros(n, m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (distm(j, i) <= bw(i))
				wtm(j, i) = pow(1 - pow(distm(j, i), 3) / pow(bw(i), 3), 3);
		}
	}
	return wtm;
}
//exponential kernel weight
// [[Rcpp::export]]
vec exp_wt_vec(vec distv, double bw)
{
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for (int i = 0; i < n; i++)
	{
		wtv(i) = exp(-distv(i) / bw);
	}
	return wtv;
}
// [[Rcpp::export]]
mat exp_wt_mat(mat distm, vec bw)
{
	int m = distm.n_cols;
	int n = distm.n_rows;
	mat wtm;
	wtm.zeros(n, m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			wtm(j, i) = exp(-distm(j, i) / bw(i));
		}
	}
	return wtm;
}
//GWR clalibration
// [[Rcpp::export]]
List gw_reg(mat x, vec y, vec w, bool hatmatrix, int focus)
{
	mat wspan(1, x.n_cols, fill::ones);
	mat xtw = trans(x % (w * wspan));
	mat xtwx = xtw * x;
	mat xtwy = trans(x) * (w % y);
	mat xtwx_inv = inv(xtwx);
	vec beta = xtwx_inv * xtwy;
	if (hatmatrix)
	{
		mat ci = xtwx_inv * xtw;
		mat s_ri = x.row(focus - 1) * ci;
		return List::create(
				Named("beta") = beta,
				Named("S_ri") = s_ri,
				Named("Ci") = ci);
	}
	else
	{
		return List::create(
				Named("beta") = beta);
	}
}
// Trace of hat matrix + trace of HH' in one function. Used in beta_se

// [[Rcpp::export]]
vec trhat2(mat S)
{
	int n_obs = S.n_rows;
	double htr = 0.0;
	double htr2 = 0.0;
	vec result(2);
	for (int i = 0; i < n_obs; i++)
	{
		htr += S(i, i);
		htr2 += sum(S.row(i) % S.row(i));
	}
	result(0) = htr;
	result(1) = htr2;
	return result;
}
// Fited values

// [[Rcpp::export]]
vec fitted(mat X, mat beta)
{
	vec fitted = sum(beta % X, 1);
	return fitted;
}

// Residuals

// [[Rcpp::export]]
vec ehat(vec y, mat X, mat beta)
{
	vec fitted = sum(beta % X, 1);
	return y - fitted;
}

// Residual sum-of squares
// [[Rcpp::export]]
double rss(vec y, mat X, mat beta)
{
	vec r = ehat(y, X, beta);
	return sum(r % r);
}

// [[Rcpp::export]]
vec gwr_diag(vec y,mat x, mat beta, mat S) {
  double ss = rss(y,x,beta);
  vec s_hat = trhat2(S);
  int n = S.n_rows;
  vec result(9);
  double AIC = n*log(ss/n)+n*log(2*datum::pi)+n+s_hat(0); //AIC
  double AICc = n*log(ss/n)+n*log(2*datum::pi)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc
  double edf = n-2*s_hat(0) + s_hat(1); //edf
  double enp = 2*s_hat(0) - s_hat(1); // enp
  double yss = sum(pow(y-mean(y),2)); //yss.g
  double r2 = 1 - ss/yss;
  double r2_adj = 1-(1-r2)*(n-1)/(edf-1);
  result(0) = AIC;
  result(1) = AICc;
  result(2) = edf;
  result(3) = enp;
  result(4) = ss;
  result(5) = r2;
  result(6) = r2_adj;
  result(7) = s_hat(0);
  result(8) = s_hat(1);
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// [[Rcpp::export]]
vec gwr_diag1(vec y, mat x, mat beta, vec s_hat)
{
	double ss = rss(y, x, beta);
	// vec s_hat = trhat2(S);
	int n = x.n_rows;
	// vec result(9);
	vec result(7);
	double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);																//AIC
	double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
	double edf = n - 2 * s_hat(0) + s_hat(1);																														//edf
	double enp = 2 * s_hat(0) - s_hat(1);																																// enp
	double yss = sum(pow(y - mean(y), 2));																															//yss.g
	double r2 = 1 - ss / yss;
	double r2_adj = 1 - (1 - r2) * (n - 1) / (edf - 1);
	result(0) = AIC;
	result(1) = AICc;
	result(2) = edf;
	result(3) = enp;
	result(4) = ss;
	result(5) = r2;
	result(6) = r2_adj;
	// result(7) = s_hat(0);
	// result(8) = s_hat(1);
	return result;
	//return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc - nice tool to give to 'optimise' to find 'best' bandwidth
// [[Rcpp::export]]
double AICc(vec y, mat x, mat beta, mat S)
{
	double ss = rss(y, x, beta);
	vec s_hat = trhat2(S);
	int n = S.n_rows;
	double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
	return AICc;
	//return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// [[Rcpp::export]]
double AICc1(vec y, mat x, mat beta, vec s_hat)
{
  double ss = rss(y, x, beta);
  // vec s_hat = trhat2(S);
  int n = x.n_rows;
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  return AICc;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc and RSS , used for function model.selection
// [[Rcpp::export]]
vec AICc_rss(vec y, mat x, mat beta, mat S)
{
	vec result(3);
	double ss = rss(y, x, beta);
	result[0] = ss;
	vec s_hat = trhat2(S);
	int n = S.n_rows;
	double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);
	double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
	result[1] = AIC;
	result[2] = AICc;
	return result;
	//return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// [[Rcpp::export]]
vec AICc_rss1(vec y, mat x, mat beta, vec s_hat)
{
  vec result(3);
  double ss = rss(y, x, beta);
  result[0] = ss;
  int n = x.n_rows;
  double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  result[1] = AIC;
  result[2] = AICc;
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

//Caculate the i row of
// [[Rcpp::export]]
mat Ci_mat(mat x, vec w)
{
	return inv(trans(x) * diagmat(w) * x) * trans(x) * diagmat(w);
}

//Scalable GWR C++ functions
// [[Rcpp::export]]
List scgwr_pre(mat x, vec y, int bw, int poly, double b0, mat g0, mat neighbour) {
  int n = x.n_rows;
  int k = x.n_cols;
  mat g0s(g0.n_cols + 1, g0.n_rows, fill::ones);
  mat g0t = trans(g0);
  for (int i = 0; i < bw; i++) {
    g0s.row(i + 1) = g0t.row(i);
  }
  g0s = trans(g0s);
  mat Mx0((poly + 1)*k*k, n, fill::zeros);
  mat My0((poly + 1)*k, n, fill::zeros);
  mat spanXnei(1, poly + 1, fill::ones);
  mat spanXtG(1, k, fill::ones);
  for (int i = 0; i < n; i++) {
    mat g(poly + 1, bw + 1, fill::ones);
    for (int p = 0; p < poly; p++) {
      g.row(p + 1) = pow(g0s.row(i), pow(2.0, poly/2.0)/pow(2.0, p + 1));
    }
    g = trans(g);
    g = g.rows(1, bw);
    mat xnei(bw, k, fill::zeros);
    vec ynei(bw, fill::zeros);
    for (int j = 0; j < bw; j++) {
      int inei = int(neighbour(i, j) - 1);
      xnei.row(j) = x.row(inei);
      ynei.row(j) = y(inei);
    }
    for (int k1 = 0; k1 < k; k1++) {
      mat XtG = xnei.col(k1) * spanXnei % g;
      for (int p = 0; p < (poly + 1); p++) {
        mat XtGX = XtG.col(p) * spanXtG % xnei;
        for (int k2 = 0; k2 < k; k2++) {
          int xindex = (k1 * (poly + 1) + p) * k + k2;
          Mx0(xindex, i) = sum(XtGX.col(k2));
        }
        int yindex = p * k + k1;
        vec XtGY = XtG.col(p) % ynei;
        My0(yindex, i) = sum(XtGY);
      }
    }
  }
  return List::create(
    Named("Mx0") = Mx0,
    Named("My0") = My0
  );
}


// [[Rcpp::export]]
double scgwr_loocv(vec target, mat x, vec y, int bw, int poly, mat Mx0, mat My0, mat XtX, mat XtY) {
  int n = x.n_rows, k = x.n_cols, poly1 = poly + 1;
  double b = target(0) * target(0), a = target(1) * target(1);
  vec R0 = vec(poly1) * b;
  for (int p = 1; p < poly1; p++) {
    R0(p) = pow(b, p + 1);
  }
  vec Rx(k*k*poly1, fill::zeros), Ry(k*poly1, fill::zeros);
  for (int p = 0; p < poly1; p++) {
    for (int k2 = 0; k2 < k; k2++) {
      for (int k1 = 0; k1 < k; k1++) {
        int xindex = k1*poly1*k + p*k + k2;
        Rx(xindex) = R0(p);
      }
      int yindex = p*k + k2;
      Ry(yindex) = R0(p);
    }
  }
  mat Mx = Rx * mat(1, n, fill::ones) % Mx0, My = Ry * mat(1, n, fill::ones) % My0;
  vec yhat(n, 1, fill::zeros);
  for (int i = 0; i < n; i++) {
    mat sumMx(k, k, fill::zeros);
    vec sumMy(k, fill::zeros);
    for (int k2 = 0; k2 < k; k2++) {
      for (int p = 0; p < poly1; p++) {
        for (int k1 = 0; k1 < k; k1++) {
          int xindex = k1*poly1*k + p*k + k2;
          sumMx(k1, k2) += Mx(xindex, i);
        }
        int yindex = p*k + k2;
        sumMy(k2) += My(yindex, i);
      }
    }
    sumMx = sumMx + a * XtX;
    sumMy = sumMy + a * XtY;
    if (det(sumMx) < 1e-10) {
      return 1e6;
    } else {
      mat beta = solve(sumMx, sumMy);
      yhat.row(i) = x.row(i) * beta;
    }
  }
  return sum((y - yhat) % (y - yhat));
}


// [[Rcpp::export]]
List scgwr_reg(mat x, vec y, int bw, int poly, mat G0, mat Mx0, mat My0, mat XtX, mat XtY, mat neighbour, vec parameters) {
  int n = x.n_rows, k = x.n_cols, poly1 = poly + 1;
  double b = parameters(0), a = parameters(1);
  /**
   * Calculate Rx, Ry, and R0.
   */
  // printf("Calculate Rx, Ry, and R0 ");
  vec R0 = vec(poly1, fill::ones) * b;
  for (int p = 1; p < poly1; p++) {
    R0(p) = pow(b, p + 1);
  }
  vec Rx(k*k*poly1, fill::zeros), Ry(k*poly1, fill::zeros);
  for (int p = 0; p < poly1; p++) {
    for (int k2 = 0; k2 < k; k2++) {
      for (int k1 = 0; k1 < k; k1++) {
        Rx(k1*poly1*k + p*k + k2) = R0(p);
      }
      Ry(p*k + k2) = R0(p);
    }
  }
  /**
  * Create G0.
  */
  // printf("Create G0 ");
  mat G0s(G0.n_cols + 1, G0.n_rows, fill::ones);
  mat G0t = trans(G0);
  G0s.rows(1, bw) = G0t.rows(0, bw - 1);
  G0s = trans(G0s);
  mat G2(n, bw + 1, fill::zeros);
  for (int p = 0; p < poly; p++) {
    G2 += pow(G0s, pow(2.0, poly/2.0)/pow(2.0, p + 1)) * R0(poly - 1);
  }
  /**
   * Calculate Mx, My.
   */
  // printf("Calculate Mx, My ");
  // mat Mx00(Mx0), My00(My0);
  for (int i = 0; i < n; i++) {
    for (int k1 = 0; k1 < k; k1++) {
      for (int p = 0; p < poly1; p++) {
        for (int k2 = 0; k2 < k; k2++) {
          Mx0((k1 * (poly + 1) + p) * k + k2, i) += x(i, k1) * x(i, k2);
        }
        My0(p * k + k1, i) += x(i, k1) * y(i);
      }
    }
  }
  mat Mx = (Rx * mat(1, n, fill::ones)) % Mx0, My = (Ry * mat(1, n, fill::ones)) % My0;
  /**
   * Regression.
   */
  // printf("Regression ");
  mat Xp(bw + 1, k * poly1, fill::zeros);
  mat rowsumG(poly1, 1, fill::ones);
  mat colsumXp(1, bw + 1, fill::zeros);
  mat spanG(1, k, fill::ones);
  mat spanX(1, poly1, fill::ones);
  mat betas(n, k, fill::zeros);
  mat betasSE(n, k, fill::zeros);
  double trS = 0.0, trStS = 0.0;
  for (int i = 0; i < n; i++) {
    /**
     * Calculate G.
     */
    mat G = mat(poly1, bw + 1, fill::ones) * R0(0);
    for (int p = 0; p < poly; p++) {
      G.row(p + 1) = pow(G0s.row(i), pow(2.0, poly/2.0)/pow(2.0, p + 1)) * R0(p);
    }
    G = trans(G);
    mat g = G * rowsumG;
    /**
     * Calculate Xp.
     */
    mat xnei(bw + 1, k, fill::zeros);
    vec ynei(bw + 1, fill::zeros);
    xnei.row(0) = x.row(i);
    ynei.row(0) = y.row(i);
    for (int j = 0; j < bw; j++) {
      int inei = int(neighbour(i, j) - 1);
      xnei.row(j+1) = x.row(inei);
      ynei.row(j+1) = y(inei);
    }
    /**
     * Calculate sumMx, sumMy.
     */
    mat sumMx(k, k, fill::zeros);
    vec sumMy(k, fill::zeros);
    for (int k2 = 0; k2 < k; k2++) {
      for (int p = 0; p < poly1; p++) {
        for (int k1 = 0; k1 < k; k1++) {
          int xindex = k1*poly1*k + p*k + k2;
          sumMx(k1, k2) += Mx(xindex, i);
        }
        int yindex = p*k + k2;
        sumMy(k2) += My(yindex, i);
      }
    }
    sumMx += a * XtX;
    sumMy += a * XtY;
    mat invMx = inv(sumMx);
    betas.row(i) = trans(invMx * sumMy);
    /**
     * Calculate Diagnostic statistics, trS and trStS.
     */
    mat StS = invMx * trans(x.row(i));
    trS += det((x.row(i) * g(0, 0)) * StS);
    mat XG2X(k, k, fill::zeros);
    for (int k1 = 0; k1 < k; k1++) {
      for (int k2 = 0; k2 < k; k2++) {
        mat Gi = G2.row(i);
        XG2X(k1, k2) = sum(xnei.col(k1) % trans(Gi % Gi + 2 * a * Gi) % xnei.col(k2)) + a * a * XtX(k1, k2);
      }
    }
    mat XX = invMx * XG2X * invMx;
    mat xi = x.row(i);
    trStS += det(sum(xi * XX * trans(xi)));
    betasSE.row(i) = trans(sqrt(XX.diag()));
  }
  return List::create(
    Named("betas") = betas,
    Named("tr.S") = trS,
    Named("tr.StS") = trStS,
    Named("betas.SE") = betasSE
  );
}

double sp_gcdist(double lon1, double lon2, double lat1, double lat2) {
  
  double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
  double w, R, a, f, D, H1, H2;
  double lat1R, lat2R, lon1R, lon2R, DE2RA;
  
  DE2RA = M_PI/180;
  a = 6378.137;              /* WGS-84 equatorial radius in km */
    f = 1.0/298.257223563;     /* WGS-84 ellipsoid flattening factor */
    
    if (fabs(lat1 - lat2) < DOUBLE_EPS) {
      if (fabs(lon1 - lon2) < DOUBLE_EPS) {
        return 0.0;
        /* Wouter Buytaert bug caught 100211 */
      } else if (fabs((fabs(lon1) + fabs(lon2)) - 360.0) < DOUBLE_EPS) {
        return 0.0;
      }
    }
    lat1R = lat1*DE2RA;
    lat2R = lat2*DE2RA;
    lon1R = lon1*DE2RA;
    lon2R = lon2*DE2RA;
    
    F = ( lat1R + lat2R )/2.0;
    G = ( lat1R - lat2R )/2.0;
    L = ( lon1R - lon2R )/2.0;
    
    /*
    printf("%g %g %g %g; %g %g %g\n",  *lon1, *lon2, *lat1, *lat2, F, G, L);
    */
    
    sinG2 = POWDI( sin( G ), 2 );
    cosG2 = POWDI( cos( G ), 2 );
    sinF2 = POWDI( sin( F ), 2 );
    cosF2 = POWDI( cos( F ), 2 );
    sinL2 = POWDI( sin( L ), 2 );
    cosL2 = POWDI( cos( L ), 2 );
    
    S = sinG2*cosL2 + cosF2*sinL2;
    C = cosG2*cosL2 + sinF2*sinL2;
    
    w = atan( sqrt( S/C ) );
    R = sqrt( S*C )/w;
    
    D = 2*w*a;
    H1 = ( 3*R - 1 )/( 2*C );
    H2 = ( 3*R + 1 )/( 2*S );
    
    return D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 ); 
}

vec sp_dists(mat dp, vec loc) {
  int N = dp.n_rows, j;
  vec dists(N, fill::zeros);
  double uout = loc(0), vout = loc(1);
  
  for (j = 0; j < N; j++) {
    dists(j) = sp_gcdist(dp(j, 0), uout, dp(j, 1), vout);
  }
  return dists;
}

// [[Rcpp::export]]
mat gw_dist(mat dp, mat rp, int focus, double p, double theta, bool longlat, bool rp_given) {
  int ndp = dp.n_rows, nrp = rp.n_rows;
  int isFocus = focus > -1;
  mat dists;
  if (p != 2 && theta != 0 && !longlat) {
    dp = coordinate_rotate(dp, theta);
    rp = coordinate_rotate(rp, theta);
  }
  if (isFocus) {
    mat prp = trans(rp.row(focus));
    if (longlat) {
      return sp_dists(dp, prp);
    } else {
      if (p == 2.0)
        return eu_dist_vec(dp, prp);
      else if(p == 1.0)
        return cd_dist_vec(dp, prp);
      else if(p == -1.0)
        return md_dist_vec(dp, prp);
      else
        return mk_dist_vec(dp, prp, p);
    }
  } else {
    if (longlat) {
      mat dists(ndp, nrp, fill::zeros);
      for (int i = 0; i < nrp; i++) {
        dists.col(i) = sp_dists(dp, trans(rp.row(i)));
      }
      return trans(dists);
    } else {
      if (p == 2.0)
        return rp_given ? eu_dist_mat(dp, rp) : eu_dist_smat(dp);
      else if (p == 1.0)
        return rp_given ? cd_dist_mat(dp, rp) : cd_dist_smat(dp);
      else if (p == -1.0)
        return rp_given ? md_dist_mat(dp, rp) : md_dist_smat(dp);
      else
        return rp_given ? mk_dist_mat(dp, rp, p) : mk_dist_smat(dp, p);
    }
  }
}

double gw_weight_gaussian(double dist, double bw) {
  return exp(pow(dist, 2)/((-2)*pow(bw, 2)));
}

double gw_weight_exponential(double dist, double bw) {
  return exp(-dist/bw);
}

double gw_weight_bisquare(double dist, double bw) {
  return dist > bw ? 0 : pow(1 - pow(dist, 2)/pow(bw, 2), 2);
}

double gw_weight_tricube(double dist, double bw) {
  return dist > bw ? 0 : pow(1 - pow(dist, 3)/pow(bw, 3), 3);
}

double gw_weight_boxcar(double dist, double bw) {
  return dist > bw ? 0 : 1;
}

typedef double (*KERNEL)(double, double);
const KERNEL GWRKernel[5] = {
  gw_weight_gaussian,
  gw_weight_exponential,
  gw_weight_bisquare,
  gw_weight_tricube,
  gw_weight_boxcar
};

// [[Rcpp::export]]
mat gw_weight(mat dist, double bw, int kernel, bool adaptive) {
  const KERNEL *kerf = GWRKernel + kernel;
  int nr = dist.n_rows, nc = dist.n_cols;
  mat w(nr, nc, fill::zeros);
  if (adaptive) {
    for (int c = 0; c < nc; c++) {
      double dn = bw / nr, fixbw = 0;
      if (dn <= 1) {
        vec vdist = sort(dist.col(c));
        fixbw = vdist(int(bw) - 1);
      } else {
        fixbw = dn * max(dist.col(c));
      }
      for (int r = 0; r < nr; r++) {
        w(r, c) = (*kerf)(dist(r, c), fixbw);
      }
    }
  } else {
    for (int c = 0; c < nc; c++) {
      for (int r = 0; r < nr; r++) {
        w(r, c) = (*kerf)(dist(r, c), bw);
      }
    }
  }
  return w;
}

// [[Rcpp::export]]
List gw_reg_all(mat x, vec y, mat dp, bool rp_given, mat rp, bool dm_given, mat dmat, bool hatmatrix, 
                double p, double theta, bool longlat, 
                double bw, int kernel, bool adaptive,
                int ngroup, int igroup) {
  int n = rp.n_rows, k = x.n_cols;
  mat betas(n, k, fill::zeros);
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  if (hatmatrix) {
    mat betasSE(n, k, fill::zeros);
    mat s_hat(1, 2, fill::zeros);
    mat qdiag(1, n, fill::zeros);
    mat rowsumSE(n, 1, fill::ones);
    // clock_t clock0 = clock(), clock1;
    for (int i = iStart; i < iEnd; i++) {
      mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
      mat w = gw_weight(d, bw, kernel, adaptive);
      mat ws(1, k, fill::ones);
      mat xtw = trans(x %(w * ws));
      mat xtwx = xtw * x;
      mat xtwy = trans(x) * (w % y);
      mat xtwx_inv = inv(xtwx);
      betas.row(i) = trans(xtwx_inv * xtwy);
      // hatmatrix
      mat ci = xtwx_inv * xtw;
      mat si = x.row(i) * ci;
      betasSE.row(i) = trans((ci % ci) * rowsumSE);
      s_hat(0) += si(0, i);
      s_hat(1) += det(si * trans(si));
      mat onei(1, n, fill::zeros);
      onei(i) = 1;
      mat p = onei - si;
      qdiag += p % p;
    }
    return List::create(
      Named("betas") = betas,
      Named("betas.SE") = betasSE,
      Named("s_hat") = s_hat,
      Named("q.diag") = qdiag
    );
  } else {
    for (int i = iStart; i < iEnd; i++) {
      mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
      mat w = gw_weight(d, bw, kernel, adaptive);
      mat ws(1, x.n_cols, fill::ones);
      mat xtw = trans(x %(w * ws));
      mat xtwx = xtw * x;
      mat xtwy = trans(x) * (w % y);
      mat xtwx_inv = inv(xtwx);
      betas.row(i) = trans(xtwx_inv * xtwy);
    }
    return List::create(
      Named("betas") = betas
    );
  }
}

// [[Rcpp::export]]
#ifdef _OPENMP
List gw_reg_all_omp(mat x, vec y, mat dp, bool rp_given, mat rp, bool dm_given, mat dmat, bool hatmatrix, 
                    double p, double theta, bool longlat, 
                    double bw, int kernel, bool adaptive,
                    int threads, int ngroup, int igroup) {
  int n = rp.n_rows, k = x.n_cols;
  mat betas(n, k, fill::zeros);
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  if (hatmatrix) {
    mat betasSE(n, k, fill::zeros);
    mat s_hat(1, 2, fill::zeros);
    mat qdiag(1, n, fill::zeros);
    mat rowsumSE(n, 1, fill::ones);
    vec s_hat1(n, fill::zeros), s_hat2(n, fill::zeros);
    int thread_nums = threads > 0 ? threads : omp_get_num_procs() - 1;
    mat qdiag_all(thread_nums, n, fill::zeros);
    bool flag_error = false;
#pragma omp parallel for num_threads(thread_nums)
    for (int i = iStart; i < iEnd; i++) {
      if (!flag_error) {
        int thread_id = omp_get_thread_num();
        mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
        mat w = gw_weight(d, bw, kernel, adaptive);
        mat ws(1, k, fill::ones);
        mat xtw = trans(x %(w * ws));
        mat xtwx = xtw * x;
        mat xtwy = trans(x) * (w % y);
        try {
          mat xtwx_inv = inv(xtwx);
          betas.row(i) = trans(xtwx_inv * xtwy);
          // hatmatrix
          mat ci = xtwx_inv * xtw;
          mat si = x.row(i) * ci;
          betasSE.row(i) = trans((ci % ci) * rowsumSE);
          // s_hat(0) += si(0, i);
          // s_hat(1) += det(si * trans(si));
          s_hat1(i) = si(0, i);
          s_hat2(i) = det(si * trans(si));
          mat onei(1, n, fill::zeros);
          onei(i) = 1;
          mat p = onei - si;
          qdiag_all.row(thread_id) += p % p;
        } catch (...) {
          flag_error = true;
        }
      }
    }
    if (flag_error) {
      throw exception("Matrix seems singular");
    } else {
      s_hat(0) = sum(s_hat1);
      s_hat(1) = sum(s_hat2);
      qdiag = mat(1, thread_nums, fill::ones) * qdiag_all;
      return List::create(
        Named("betas") = betas,
        Named("betas.SE") = betasSE,
        Named("s_hat") = s_hat,
        Named("q.diag") = qdiag
      );
    }
  } else {
    bool flag_error = false;
    for (int i = iStart; i < iEnd; i++) {
      if (!flag_error) {
        mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
        mat w = gw_weight(d, bw, kernel, adaptive);
        mat ws(1, x.n_cols, fill::ones);
        mat xtw = trans(x %(w * ws));
        mat xtwx = xtw * x;
        mat xtwy = trans(x) * (w % y);
        try {
          mat xtwx_inv = inv(xtwx);
          betas.row(i) = trans(xtwx_inv * xtwy);
        } catch (...) {
          flag_error = true;
        }
      }
    }
    if (flag_error) {
      throw exception("Matrix seems singular.");
    } else {
      return List::create(
        Named("betas") = betas
      );
    }
  }
}
#endif

// [[Rcpp::export]]
double gw_cv_all(mat x, vec y, mat dp, bool dm_given, mat dmat, 
                 double p, double theta, bool longlat, 
                 double bw, int kernel, bool adaptive,
                 int ngroup, int igroup) {
  int n = dp.n_rows;
  double cv = 0.0;
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  for (int i = iStart; i < iEnd; i++) {
    mat d = dm_given ? dmat.col(i) : gw_dist(dp, dp, i, p, theta, longlat, false);
    mat w = gw_weight(d, bw, kernel, adaptive);
    w(i, 0) = 0.0;
    mat ws(1, x.n_cols, fill::ones);
    mat xtw = trans(x %(w * ws));
    mat xtwx = xtw * x;
    mat xtwy = trans(x) * (w % y);
    mat xtwx_inv = inv(xtwx);
    mat betas = xtwx_inv * xtwy;
    double res = y(i) - det(x.row(i) * betas);
    cv += res * res;
  }
  return cv;
}

// [[Rcpp::export]]
#ifdef _OPENMP
double gw_cv_all_omp(mat x, vec y, mat dp, bool dm_given, mat dmat, 
                     double p, double theta, bool longlat, 
                     double bw, int kernel, bool adaptive,
                     int threads, int ngroup, int igroup) {
  int n = dp.n_rows;
  int thread_nums = omp_get_num_procs() - 1;
  vec cv(thread_nums, fill::zeros);
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  bool flag_error = false;
#pragma omp parallel for num_threads(thread_nums)
  for (int i = iStart; i < iEnd; i++) {
    if (!flag_error) {
      int thread_id = threads > 0 ? threads : omp_get_thread_num();
      mat d = dm_given ? dmat.col(i) : gw_dist(dp, dp, i, p, theta, longlat, false);
      mat w = gw_weight(d, bw, kernel, adaptive);
      w(i, 0) = 0.0;
      mat ws(1, x.n_cols, fill::ones);
      mat xtw = trans(x %(w * ws));
      mat xtwx = xtw * x;
      mat xtwy = trans(x) * (w % y);
      try {
        mat xtwx_inv = inv(xtwx);
        mat betas = xtwx_inv * xtwy;
        double res = y(i) - det(x.row(i) * betas);
        cv(thread_id) += res * res;
      } catch (...) {
        flag_error = true;
      }
    }
  }
  if (flag_error) {
    throw exception("Matrix seems singular.");
  }
  return sum(cv);
}
#endif

// [[Rcpp::export]]
vec gw_local_r2(mat dp, vec dybar2, vec dyhat2, bool dm_given, mat dmat, double p, double theta, bool longlat, double bw, int kernel, bool adaptive) {
  int n = dp.n_rows;
  vec localR2(n, fill::zeros);
  for (int i = 0; i < n; i++) {
    mat d = dm_given ? dmat.col(i) : gw_dist(dp, dp, i, p, theta, longlat, false);
    mat w = gw_weight(d, bw, kernel, adaptive);
    double tss = sum(dybar2 % w);
    double rss = sum(dyhat2 % w);
    localR2(i) = (tss - rss) / tss;
  }
  return localR2;
}

//gwr.mixed fast code
// FE EDITS BELOW BY FIONA.H.EVANS@GMAIL.COM
// [[Rcpp::export]]
double BIC(vec y, mat x, mat beta, vec s_hat)
{
  double ss = rss(y, x, beta);
  double n = (double)x.n_rows;
  double BIC = n * log(ss / n) + n * log(2 * datum::pi) + log(n) * s_hat(0);
  return BIC;
}
// For debugging
void printVec(vec v) {
  int n = v.size();
  n = 10;
  
  for (int i=0; i<n; i++) {
    Rprintf("%f ", v(i));
  }
  Rprintf("\n");
}
void printMat(mat m) {
  uword n = m.n_rows;
  uword p = m.n_cols;
  
  n = 10;
  if (m.n_rows < n) 
  {
     n = m.n_rows;
  } 
  for (uword i=0; i<n; i++) {
    for (uword j=0; j<p; j++) {
      Rprintf("%f ", m(i, j));
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}
// Boxcar weights 
// [[Rcpp::export]]
vec box_wt_vec(vec distv, double bw)
{
  int n = distv.n_elem;
  vec wtv(n, fill::zeros);
  
  uvec u = find(distv <= bw);
  wtv.elem(u).fill(1);
  
  return wtv;
}
// Boxcar adaptive weights
// [[Rcpp::export]]
vec box_wt_vec_ad(vec distv, double bw)
{
  int n = distv.n_elem;
  vec wtv;
  wtv.zeros(n);
  double bwd = bw;
  if (bw >= n) bwd = n;
  
  // equivalent to R function rank(distv, ties.method='first')
  uvec rnk1 = sort_index(distv) + 1;
  uvec rnk = sort_index(rnk1) + 1;
  
  uvec u = find(rnk <= bwd);
  wtv.elem(u).fill(1);
  
  return wtv;
}
// Gaussian adaptive weights
// [[Rcpp::export]]
vec gau_wt_vec_ad(vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    uvec rnk1 = sort_index(distv) + 1;
    uvec rnk = sort_index(rnk1) + 1;
    
    uvec u = find(rnk == bw);
    bwd = distv(u(0));
  }
  
  vec wtv = exp(pow(distv, 2) / ((-2) * pow(bwd, 2)));
  
  return wtv;
}
// Bisquare adaptive weights
// [[Rcpp::export]]
vec bis_wt_vec_ad(vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    uvec rnk1 = sort_index(distv) + 1;
    uvec rnk = sort_index(rnk1) + 1;
    
    uvec u = find(rnk == bw);
    bwd = distv(u(0));
  }
  
  vec wtv = bisq_wt_vec(distv, bwd);
  return wtv;
}
// Tricube adaptive weight
// [[Rcpp::export]]
vec tri_wt_vec_ad(vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    uvec rnk1 = sort_index(distv) + 1;
    uvec rnk = sort_index(rnk1) + 1;
    
    uvec u = find(rnk == bw);
    bwd = distv(u(0));
  }
  
  vec wtv = tri_wt_vec(distv, bwd);
  
  return wtv;
}
// Exponential adaptive weights
// [[Rcpp::export]]
vec exp_wt_vec_ad(vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    uvec rnk1 = sort_index(distv) + 1;
    uvec rnk = sort_index(rnk1) + 1;
    
    uvec u = find(rnk == bw);
    bwd = distv(u(0));
  }
  
  vec wtv = exp_wt_vec(distv, bwd);
  
  return wtv;
}
// For use in gw_weight
enum string_code{ga, bi, tr, bo, ex};
string_code hashit (std::string const& inString) {
  if (inString == "gaussian")
  {
    return ga;
  }
  if (inString == "bisquare")
  {
    return bi;
  }
  if (inString == "tricube")
  {
    return tr;
  }
  if (inString == "boxcar")
  {
    return bo;
  }
  if (inString == "exponential") 
  {
    return ex;
  }
  return ga;
}
// Geographic weights (gw.weight equivalent)
// [[Rcpp::export]]
vec gw_weight(vec vdist, double bw, std::string kernel, bool adaptive)
{
  vec wv;
  if (adaptive) switch(hashit(kernel)){
  case ga: wv = gau_wt_vec_ad (vdist, bw);
  case bi: wv = bis_wt_vec_ad (vdist, bw);
  case tr: wv = tri_wt_vec_ad (vdist, bw);
  case bo: wv = box_wt_vec_ad (vdist, bw); 
  case ex: wv = exp_wt_vec_ad (vdist, bw);
  }
  else switch(hashit(kernel)){
  case ga: wv = gauss_wt_vec (vdist, bw);
  case bi: wv = bisq_wt_vec (vdist, bw);
  case tr: wv = tri_wt_vec (vdist, bw);
  case bo: wv = box_wt_vec (vdist, bw);
  case ex: wv = exp_wt_vec (vdist, bw);
  }
  return wv;
}
// GWR calibration, returns betas only
// [[Rcpp::export]]
vec gw_reg_2(mat x, vec y, vec w)
{
  mat wspan(1, x.n_cols, fill::ones);
  mat xtw = trans(x % (w * wspan));
  mat xtwx = xtw * x;
  mat xtwy = trans(x) * (w % y);
  mat xtwx_inv = inv(xtwx);
  vec beta = xtwx_inv * xtwy;
  
  return beta;
}
// C++ version of gwr.q, used in gwr.mixed
// [[Rcpp::export]]
mat gwr_q(mat x, vec y, 
                mat dMat, double bw, std::string kernel, bool adaptive)
{
  // loc and out.loc only used to create distances
  int n =  dMat.n_cols;  // loc.n_rows() 
  int m =  x.n_cols;
  mat beta(n, m);
  vec distv;
  vec w;
  
  for (int i = 0; i < n; i++) {
    distv = dMat.col(i);
    w = gw_weight(distv, bw, kernel, adaptive);
    beta.row(i) = gw_reg_2(x, y, w);
  }
  
  return beta;
}
// [[Rcpp::export]]
vec e_vec(int m, int n){
  vec e = linspace(0, n-1, n);
  vec ret(n, fill::zeros);
  uvec u = find(e == m);
  ret.elem(u).fill(1);
  return ret;
}
// [[Rcpp::export]]
double gwr_mixed_trace(mat x1, mat x2, vec y, 
                       mat dMat, double bw, std::string kernel, bool adaptive){
  int i;
  int n = x1.n_rows;
  int nc2 = x2.n_cols;
  mat mtemp, model1, model2;
  mat x3(n, nc2);
  vec y2(n);
  vec y3;
  vec hii(n, fill::zeros);
  mat m(1,n);
  double s1, s2;
  
  for (i = 0; i < nc2; i++) {
    mtemp = gwr_q(x1, x2.col(i), dMat, bw, kernel, adaptive);
    x3.col(i) = x2.col(i) - fitted(x1, mtemp);
  }
  
  // The following works but is slow
  for (i = 0; i < n; i++) {
    mtemp = gwr_q(x1, e_vec(i, n), dMat, bw, kernel, adaptive); // length n x nc2
    y2 = e_vec(i, n) - fitted(x1, mtemp); // length n
    model2 = gwr_q(x3, y2, dMat, 100000, "boxcar", true);
    y3 = e_vec(i, n) - fitted(x2, model2);
    model1 = gwr_q(x1, y3, dMat.col(i), bw, kernel, adaptive); // 1 x 1 matrix
    model2 = gwr_q(x3, y2, dMat.col(i), 100000, "boxcar", true); // n x nc2
    s1 = fitted(x1.row(i), model1)(0);  // vector with one element
    s2 = fitted(x2.row(i), model2)(0);  // vector with one element
    hii(i) = s1 + s2;
  }
  return sum(hii);
}
// [[Rcpp::export]]
List gwr_mixed_2(mat x1, mat x2, vec y, 
                       mat dMat, mat dMat_rp,
                       double bw, std::string kernel, bool adaptive){
  int i;
  int n = x1.n_rows;
  int nc2 = x2.n_cols;
  mat mtemp, model1, model2;
  mat x3(n, nc2);
  vec y2(n);
  vec hii(n, fill::zeros);
  mat m(1,n);  
  for (i = 0; i < nc2; i++) {
    mtemp = gwr_q(x1, x2.col(i), dMat, bw, kernel, adaptive);
    x3.col(i) = x2.col(i) - fitted(x1, mtemp);
  }
  
  mtemp = gwr_q(x1, y, dMat, bw, kernel, adaptive);
  y2 = y - fitted(x1, mtemp);
  model2 = gwr_q(x3, y2, dMat, 100000, "boxcar", true);
  
  model1 = gwr_q(x1, y-fitted(x2, model2), dMat_rp, bw, kernel, adaptive);
  model2 = gwr_q(x3, y2, dMat_rp, 100000, "boxcar", true);
  
  return List::create(
    Named("local") = model1,
    Named("global") = model2
  );
}
