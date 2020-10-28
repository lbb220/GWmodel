#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP GWmodel_AICc(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_AICc1(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_AICc_rss(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_AICc_rss1(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_bisq_wt_mat(SEXP, SEXP);
extern SEXP GWmodel_bisq_wt_vec(SEXP, SEXP);
extern SEXP GWmodel_cd_dist_mat(SEXP, SEXP);
extern SEXP GWmodel_cd_dist_smat(SEXP);
extern SEXP GWmodel_cd_dist_vec(SEXP, SEXP);
extern SEXP GWmodel_Ci_mat(SEXP, SEXP);
extern SEXP GWmodel_coordinate_rotate(SEXP, SEXP);
extern SEXP GWmodel_ehat(SEXP, SEXP, SEXP);
extern SEXP GWmodel_eu_dist_mat(SEXP, SEXP);
extern SEXP GWmodel_eu_dist_smat(SEXP);
extern SEXP GWmodel_eu_dist_vec(SEXP, SEXP);
extern SEXP GWmodel_exp_wt_mat(SEXP, SEXP);
extern SEXP GWmodel_exp_wt_vec(SEXP, SEXP);
extern SEXP GWmodel_fitted(SEXP, SEXP);
extern SEXP GWmodel_gauss_wt_mat(SEXP, SEXP);
extern SEXP GWmodel_gauss_wt_vec(SEXP, SEXP);
extern SEXP GWmodel_gw_reg(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gwr_diag(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gwr_diag1(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_md_dist_mat(SEXP, SEXP);
extern SEXP GWmodel_md_dist_smat(SEXP);
extern SEXP GWmodel_md_dist_vec(SEXP, SEXP);
extern SEXP GWmodel_mk_dist_mat(SEXP, SEXP, SEXP);
extern SEXP GWmodel_mk_dist_smat(SEXP, SEXP);
extern SEXP GWmodel_mk_dist_vec(SEXP, SEXP, SEXP);
extern SEXP GWmodel_rss(SEXP, SEXP, SEXP);
extern SEXP GWmodel_tri_wt_mat(SEXP, SEXP);
extern SEXP GWmodel_tri_wt_vec(SEXP, SEXP);
extern SEXP GWmodel_scgwr_pre(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_scgwr_loocv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_scgwr_reg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_dist(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_reg_all(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_reg_all_omp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_cv_all(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_cv_all_omp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_local_r2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_reg_cuda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gw_cv_all_cuda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_BIC(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_box_wt_vec(SEXP, SEXP);
extern SEXP GWmodel_box_wt_vec_ad(SEXP, SEXP);
extern SEXP GWmodel_gau_wt_vec_ad(SEXP, SEXP);
extern SEXP GWmodel_bis_wt_vec_ad(SEXP, SEXP);
extern SEXP GWmodel_tri_wt_vec_ad(SEXP, SEXP);
extern SEXP GWmodel_exp_wt_vec_ad(SEXP, SEXP);
extern SEXP GWmodel_gw_reg_2(SEXP, SEXP, SEXP);
extern SEXP GWmodel_gwr_q(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_e_vec(SEXP, SEXP);
extern SEXP GWmodel_gwr_mixed_trace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gwr_mixed_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GWmodel_AICc",              (DL_FUNC) &GWmodel_AICc,              4},
    {"GWmodel_AICc1",             (DL_FUNC) &GWmodel_AICc1,             4},
    {"GWmodel_AICc_rss",          (DL_FUNC) &GWmodel_AICc_rss,          4},
    {"GWmodel_AICc_rss1",         (DL_FUNC) &GWmodel_AICc_rss1,         4},
    {"GWmodel_bisq_wt_mat",       (DL_FUNC) &GWmodel_bisq_wt_mat,       2},
    {"GWmodel_bisq_wt_vec",       (DL_FUNC) &GWmodel_bisq_wt_vec,       2},
    {"GWmodel_cd_dist_mat",       (DL_FUNC) &GWmodel_cd_dist_mat,       2},
    {"GWmodel_cd_dist_smat",      (DL_FUNC) &GWmodel_cd_dist_smat,      1},
    {"GWmodel_cd_dist_vec",       (DL_FUNC) &GWmodel_cd_dist_vec,       2},
    {"GWmodel_Ci_mat",            (DL_FUNC) &GWmodel_Ci_mat,            2},
    {"GWmodel_coordinate_rotate", (DL_FUNC) &GWmodel_coordinate_rotate, 2},
    {"GWmodel_ehat",              (DL_FUNC) &GWmodel_ehat,              3},
    {"GWmodel_eu_dist_mat",       (DL_FUNC) &GWmodel_eu_dist_mat,       2},
    {"GWmodel_eu_dist_smat",      (DL_FUNC) &GWmodel_eu_dist_smat,      1},
    {"GWmodel_eu_dist_vec",       (DL_FUNC) &GWmodel_eu_dist_vec,       2},
    {"GWmodel_exp_wt_mat",        (DL_FUNC) &GWmodel_exp_wt_mat,        2},
    {"GWmodel_exp_wt_vec",        (DL_FUNC) &GWmodel_exp_wt_vec,        2},
    {"GWmodel_fitted",            (DL_FUNC) &GWmodel_fitted,            2},
    {"GWmodel_gauss_wt_mat",      (DL_FUNC) &GWmodel_gauss_wt_mat,      2},
    {"GWmodel_gauss_wt_vec",      (DL_FUNC) &GWmodel_gauss_wt_vec,      2},
    {"GWmodel_gw_reg",            (DL_FUNC) &GWmodel_gw_reg,            5},
    {"GWmodel_gwr_diag",          (DL_FUNC) &GWmodel_gwr_diag,          4},
    {"GWmodel_gwr_diag1",          (DL_FUNC) &GWmodel_gwr_diag1,        4},
    {"GWmodel_md_dist_mat",       (DL_FUNC) &GWmodel_md_dist_mat,       2},
    {"GWmodel_md_dist_smat",      (DL_FUNC) &GWmodel_md_dist_smat,      1},
    {"GWmodel_md_dist_vec",       (DL_FUNC) &GWmodel_md_dist_vec,       2},
    {"GWmodel_mk_dist_mat",       (DL_FUNC) &GWmodel_mk_dist_mat,       3},
	{"GWmodel_mk_dist_smat",      (DL_FUNC) &GWmodel_mk_dist_smat,      2},
    {"GWmodel_mk_dist_vec",       (DL_FUNC) &GWmodel_mk_dist_vec,       3},
    {"GWmodel_rss",               (DL_FUNC) &GWmodel_rss,               3},
    {"GWmodel_tri_wt_mat",        (DL_FUNC) &GWmodel_tri_wt_mat,        2},
    {"GWmodel_tri_wt_vec",        (DL_FUNC) &GWmodel_tri_wt_vec,        2},
    {"GWmodel_scgwr_pre",         (DL_FUNC) &GWmodel_scgwr_pre,         7},
    {"GWmodel_scgwr_loocv",       (DL_FUNC) &GWmodel_scgwr_loocv,       9},
    {"GWmodel_scgwr_reg",         (DL_FUNC) &GWmodel_scgwr_reg,         11},
    {"GWmodel_gw_dist",           (DL_FUNC) &GWmodel_gw_dist,           7},
    {"GWmodel_gw_reg_all",        (DL_FUNC) &GWmodel_gw_reg_all,        16},
    {"GWmodel_gw_reg_all_omp",    (DL_FUNC) &GWmodel_gw_reg_all_omp,    17},
    {"GWmodel_gw_cv_all",         (DL_FUNC) &GWmodel_gw_cv_all,         13},
    {"GWmodel_gw_cv_all_omp",     (DL_FUNC) &GWmodel_gw_cv_all_omp,     14},
    {"GWmodel_gw_local_r2",       (DL_FUNC) &GWmodel_gw_local_r2,       11},
    {"GWmodel_gw_reg_cuda",       (DL_FUNC) &GWmodel_gw_reg_cuda,       16},
    {"GWmodel_gw_cv_all_cuda",    (DL_FUNC) &GWmodel_gw_cv_all_cuda,    13},
    {"GWmodel_BIC",               (DL_FUNC) &GWmodel_BIC,                4},
    {"GWmodel_box_wt_vec",        (DL_FUNC) &GWmodel_box_wt_vec,         2},
    {"GWmodel_box_wt_vec_ad",     (DL_FUNC) &GWmodel_box_wt_vec_ad,      2},
    {"GWmodel_gau_wt_vec_ad",     (DL_FUNC) &GWmodel_gau_wt_vec_ad,      2},
    {"GWmodel_bis_wt_vec_ad",     (DL_FUNC) &GWmodel_bis_wt_vec_ad,      2},
    {"GWmodel_tri_wt_vec_ad",     (DL_FUNC) &GWmodel_tri_wt_vec_ad,      2},
    {"GWmodel_exp_wt_vec_ad",     (DL_FUNC) &GWmodel_exp_wt_vec_ad,      2},
    {"GWmodel_gw_reg_2",          (DL_FUNC) &GWmodel_gw_reg_2,           3},
    {"GWmodel_gwr_q",             (DL_FUNC) &GWmodel_gwr_q,              6},
    {"GWmodel_e_vec",             (DL_FUNC) &GWmodel_e_vec,              2},
    {"GWmodel_gwr_mixed_trace",   (DL_FUNC) &GWmodel_gwr_mixed_trace,    7},
    {"GWmodel_gwr_mixed_2",       (DL_FUNC) &GWmodel_gwr_mixed_2,        8},
    {NULL, NULL, 0}
};

void R_init_GWmodel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
