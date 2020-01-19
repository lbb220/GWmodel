#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP GWmodel_AICc(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_AICc_rss(SEXP, SEXP, SEXP, SEXP);
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
extern SEXP GWmodel_AICc1(SEXP, SEXP, SEXP, SEXP);
extern SEXP GWmodel_gwr_diag1(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GWmodel_AICc",              (DL_FUNC) &GWmodel_AICc,              4},
    {"GWmodel_AICc_rss",          (DL_FUNC) &GWmodel_AICc_rss,          4},
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
	{"GWmodel_AICc1",              (DL_FUNC) &GWmodel_AICc1,             4},
	{"GWmodel_gwr_diag1",          (DL_FUNC) &GWmodel_gwr_diag1,        4},
    {NULL, NULL, 0}
};

void R_init_GWmodel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
