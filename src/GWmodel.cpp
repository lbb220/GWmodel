// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------


//distance matrix calculation
//coords must be a matrix with 2 columns
// [[Rcpp::export]]
mat coordinate_rotate(mat coords, double theta){
	int n = coords.n_rows;
    mat rotated_coords(n,2);
    rotated_coords.col(0) = coords.col(0)*cos(theta)-coords.col(1)*sin(theta);	
	rotated_coords.col(1) = coords.col(0)*sin(theta)+coords.col(1)*cos(theta);
    return rotated_coords;	
}

//Eudclidean distance matrix
// [[Rcpp::export]]
mat eu_dist_mat(mat in_locs, mat out_locs){
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat eu_dist(n_in,n_out);
	for(int i=0;i<n_in;i++){
		for(int j=0;j<n_out;j++){
			eu_dist(i,j) = sqrt(sum(pow(in_locs.row(i) - out_locs.row(j),2)));
		}
	}
	return eu_dist;
}
//symmetrical distance matrix
// [[Rcpp::export]]
mat eu_dist_smat(mat in_locs){
	int n = in_locs.n_rows;
	mat eu_dist(n,n);
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			eu_dist(i,j) = sqrt(sum(pow(in_locs.row(i) - in_locs.row(j),2)));
			eu_dist(j,i) = eu_dist(i,j);
		}
	}
	return eu_dist;
}
// [[Rcpp::export]]
vec eu_dist_vec(mat in_locs, vec out_loc){
	int n_in = in_locs.n_rows;
	vec eu_dist(n_in);
	for(int i=0;i<n_in;i++){
			eu_dist(i) = sqrt(sum(pow(in_locs.row(i) - trans(out_loc),2)));
	}
	return eu_dist;
}



//Manhattan distance matrix
// [[Rcpp::export]]
mat md_dist_mat(mat in_locs, mat out_locs){
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat md_dist(n_in,n_out);
	for(int i=0;i<n_in;i++){
		for(int j=0;j<n_out;j++){
			md_dist(i,j) = sum(abs(in_locs.row(i) - out_locs.row(j)));
		}
	}
	return md_dist;
}

//symmetrical distance matrix
// [[Rcpp::export]]
mat md_dist_smat(mat in_locs){
	int n = in_locs.n_rows;
	mat md_dist(n,n);
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			md_dist(i,j) = sum(abs(in_locs.row(i) - in_locs.row(j)));
			md_dist(j,i) = md_dist(i,j);
		}
	}
	return md_dist;
}
// [[Rcpp::export]]
vec md_dist_vec(mat in_locs, vec out_loc){
	int n_in = in_locs.n_rows;
	vec md_dist(n_in);
	for(int i=0;i<n_in;i++){
		md_dist(i) = sum(abs(in_locs.row(i) - trans(out_loc)));
	}
	return md_dist;
}

//Chebyshev distance matrix
// [[Rcpp::export]]
mat cd_dist_mat(mat in_locs, mat out_locs){
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat cd_dist(n_in,n_out);
	for(int i=0;i<n_in;i++){
		for(int j=i;j<n_out;j++){
			cd_dist(i,j) = max(abs(in_locs.row(i) - out_locs.row(j)));
			cd_dist(j,i) = cd_dist(i,j);
		}
	}
	return cd_dist;
}

//symmetrical distance matrix
// [[Rcpp::export]]
mat cd_dist_smat(mat in_locs){
	int n = in_locs.n_rows;
	mat cd_dist(n,n);
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			cd_dist(i,j) = max(abs(in_locs.row(i) - in_locs.row(j)));
			cd_dist(j,i) = cd_dist(i,j);
		}
	}
	return cd_dist;
}
// [[Rcpp::export]]
vec cd_dist_vec(mat in_locs, vec out_loc){
	int n_in = in_locs.n_rows;
	vec cd_dist(n_in);
	for(int i=0;i<n_in;i++){
		cd_dist(i) = max(abs(in_locs.row(i) - trans(out_loc)));
	}
	return cd_dist;
}

//Minkowski distance matrix
// [[Rcpp::export]]
mat mk_dist_mat(mat in_locs, mat out_locs,double p){
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat mk_dist(n_in,n_out);
	for(int i=0;i<n_in;i++){
		for(int j=0;j<n_out;j++){
			mk_dist(i,j) = pow(sum(pow(abs(in_locs.row(i) - out_locs.row(j)),p)),1.0/p);
		}
	}
	return mk_dist;
}
//sqrt(sum(pow(in_locs.row(i) - trans(out_loc),2)))
//symmetrical distance matrix
// [[Rcpp::export]]
mat mk_dist_smat(mat in_locs,double p){
	int n = in_locs.n_rows;
	mat mk_dist(n,n);
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			mk_dist(i,j) = pow(sum(pow(abs(in_locs.row(i) - in_locs.row(j)),p)),1.0/p);
			mk_dist(j,i) = mk_dist(i,j);
		}
	}
	return mk_dist;
}

// [[Rcpp::export]]
vec mk_dist_vec(mat in_locs, vec out_loc, double p){
	int n_in = in_locs.n_rows;
	vec mk_dist(n_in);
	for(int i=0;i<n_in;i++){
		mk_dist(i) = pow(sum(pow(abs(in_locs.row(i) - trans(out_loc)),p)),1.0/p);
	}
	return mk_dist;
}
//Weight matrix
//Bisuqare weight
// [[Rcpp::export]]
vec bisq_wt_vec(vec distv, double bw){
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for(int i=0; i<n; i++){
		if(distv(i) <= bw)
			wtv(i) = pow(1 - pow(distv(i)/bw,2),2);
	}
	return wtv;
}
//Calculated by column, the length of bw must be the same the number of columns of distm
// [[Rcpp::export]]
mat bisq_wt_mat(mat distm, vec bw){
	int m = distm.n_cols;
	int n = distm.n_rows; 
	mat wtm;
	wtm.zeros(n,m);
	for(int i=0; i<m; i++){
		for(int j=0;j<n;j++){
			if(distm(j,i) <= bw(i))
				wtm(j,i) = pow(1 - pow(distm(j,i)/bw(i),2),2);
		}
		
	}
	return wtm;
}

//Gaussian weight
// [[Rcpp::export]]
vec gauss_wt_vec(vec distv, double bw){
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for(int i=0; i<n; i++){
		wtv(i) = exp(pow(distv(i),2)/((-2)*pow(bw,2)));
	}
	return wtv;
}
// [[Rcpp::export]]
mat gauss_wt_mat(mat distm, vec bw){
	int m = distm.n_cols;
	int n = distm.n_rows; 
	mat wtm;
	wtm.zeros(n,m);
	for(int i=0; i<m; i++){
		for(int j=0;j<n;j++){
			wtm(j,i) = exp(pow(distm(j,i),2)/((-2)*pow(bw(i),2)));
		}
		
	}
	return wtm;
}

//Tricube weight
// [[Rcpp::export]]
vec tri_wt_vec(vec distv, double bw){
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for(int i=0; i<n; i++){
		if(distv(i) <= bw)
			wtv(i) = pow(1 - pow(distv(i),3)/pow(bw,3),3);
	}
	return wtv;
}
// [[Rcpp::export]]
mat tri_wt_mat(mat distm, vec bw){
	int m = distm.n_cols;
	int n = distm.n_rows; 
	mat wtm;
	wtm.zeros(n,m);
	for(int i=0; i<m; i++){
		for(int j=0;j<n;j++){
			if(distm(j,i) <= bw(i))
				wtm(j,i) = pow(1 - pow(distm(j,i),3)/pow(bw(i),3),3);
		}
		
	}
	return wtm;
}
//exponential kernel weight
// [[Rcpp::export]]
vec exp_wt_vec(vec distv, double bw){
	int n = distv.n_elem;
	vec wtv;
	wtv.zeros(n);
	for(int i=0; i<n; i++){
		wtv(i) = exp(-distv(i)/bw);
	}
	return wtv;
}
// [[Rcpp::export]]
mat exp_wt_mat(mat distm, vec bw){
	int m = distm.n_cols;
	int n = distm.n_rows; 
	mat wtm;
	wtm.zeros(n,m);
	for(int i=0; i<m; i++){
		for(int j=0;j<n;j++){
			wtm(j,i) = exp(-distm(j,i)/bw(i));
		}
		
	}
	return wtm;
}
//GWR clalibration
// [[Rcpp::export]]
List gw_reg(mat x, vec y, vec w, bool hatmatrix, int focus){
	mat xtw = trans(x)*diagmat(w);
	vec beta = solve(xtw*x,xtw*y);
	if(hatmatrix){
		mat xtwx_inv = inv(xtw*x);
		mat ci = xtwx_inv*xtw;
		mat s_ri = x.row(focus-1)*ci;
		return List::create(
		Named("beta") = beta,
		Named("S_ri") = s_ri,
		Named("Ci") = ci
		);
	}
	else{
		return List::create(
		Named("beta") = beta
		);
	}
}
// Trace of hat matrix + trace of HH' in one function. Used in beta_se

// [[Rcpp::export]]
vec trhat2(mat S) {
  int n_obs = S.n_rows;
  double htr = 0.0;
  double htr2 = 0.0;
  vec result(2);
  for (int i = 0; i < n_obs; i++) {
    htr  += S(i,i);
    htr2 += sum(S.row(i) % S.row(i));
  }
  result(0) = htr;
  result(1) = htr2;
  return result;
}
// Fited values

// [[Rcpp::export]]
vec fitted(mat X, mat beta) {
  vec fitted = sum(beta % X, 1);
  return fitted;
}

// Residuals

// [[Rcpp::export]]
vec ehat(vec y, mat X, mat beta) {
  vec fitted = sum(beta % X, 1);
  return y - fitted;
}

// Residual sum-of squares
// [[Rcpp::export]]
double rss(vec y, mat X, mat beta) {
  vec r = ehat(y, X, beta);
  return sum(r%r);
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

// return the AICc - nice tool to give to 'optimise' to find 'best' bandwidth
// [[Rcpp::export]]
double AICc(vec y,mat x, mat beta, mat S) {
  double ss = rss(y,x,beta);
  vec s_hat = trhat2(S);
  int n = S.n_rows;
  double AICc = n*log(ss/n)+n*log(2*datum::pi)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc
  return AICc;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc and RSS , used for function model.selection
// [[Rcpp::export]]
vec AICc_rss(vec y,mat x, mat beta, mat S) {
  vec result(3);
  double ss = rss(y,x,beta);
  result[0] = ss;
  vec s_hat = trhat2(S);
  int n = S.n_rows;
  double AIC= n*log(ss/n)+n*log(2*datum::pi)+n+s_hat(0);
  double AICc = n*log(ss/n)+n*log(2*datum::pi)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc
  result[1] = AIC;
  result[2] = AICc;
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}


//Caculate the i row of 
// [[Rcpp::export]]
mat Ci_mat(mat x, vec w) {
	return inv(trans(x)*diagmat(w)*x)*trans(x)*diagmat(w);
}





























