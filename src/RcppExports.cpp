// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gibbs_sampler_with_corr
List gibbs_sampler_with_corr(arma::mat& Y, arma::mat& Sd_hat, arma::mat& trait_corr, int& N, arma::mat& B, double& sigma, arma::vec& sigma1, arma::vec& sigma0, arma::vec& p, arma::mat& A, arma::umat& Z, double& alpha_B, double& beta_B, arma::vec& alpha_0, arma::vec& alpha_1, arma::vec& beta_0, arma::vec& beta_1, arma::vec& a, arma::vec& b, arma::mat& Lambda3, arma::mat& Lambda_inv3);
RcppExport SEXP _FLOWMR_gibbs_sampler_with_corr(SEXP YSEXP, SEXP Sd_hatSEXP, SEXP trait_corrSEXP, SEXP NSEXP, SEXP BSEXP, SEXP sigmaSEXP, SEXP sigma1SEXP, SEXP sigma0SEXP, SEXP pSEXP, SEXP ASEXP, SEXP ZSEXP, SEXP alpha_BSEXP, SEXP beta_BSEXP, SEXP alpha_0SEXP, SEXP alpha_1SEXP, SEXP beta_0SEXP, SEXP beta_1SEXP, SEXP aSEXP, SEXP bSEXP, SEXP Lambda3SEXP, SEXP Lambda_inv3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sd_hat(Sd_hatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trait_corr(trait_corrSEXP);
    Rcpp::traits::input_parameter< int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma1(sigma1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::umat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double& >::type alpha_B(alpha_BSEXP);
    Rcpp::traits::input_parameter< double& >::type beta_B(beta_BSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type alpha_0(alpha_0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type alpha_1(alpha_1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_1(beta_1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Lambda3(Lambda3SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Lambda_inv3(Lambda_inv3SEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_sampler_with_corr(Y, Sd_hat, trait_corr, N, B, sigma, sigma1, sigma0, p, A, Z, alpha_B, beta_B, alpha_0, alpha_1, beta_0, beta_1, a, b, Lambda3, Lambda_inv3));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_sampler
List gibbs_sampler(arma::mat& Y, arma::mat& Sd_hat, int& N, arma::mat& B, double& sigma, arma::vec& sigma1, arma::vec& sigma0, arma::vec& p, arma::mat& A, arma::umat& Z, double& alpha_B, double& beta_B, arma::vec& alpha_0, arma::vec& alpha_1, arma::vec& beta_0, arma::vec& beta_1, arma::vec& a, arma::vec& b);
RcppExport SEXP _FLOWMR_gibbs_sampler(SEXP YSEXP, SEXP Sd_hatSEXP, SEXP NSEXP, SEXP BSEXP, SEXP sigmaSEXP, SEXP sigma1SEXP, SEXP sigma0SEXP, SEXP pSEXP, SEXP ASEXP, SEXP ZSEXP, SEXP alpha_BSEXP, SEXP beta_BSEXP, SEXP alpha_0SEXP, SEXP alpha_1SEXP, SEXP beta_0SEXP, SEXP beta_1SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sd_hat(Sd_hatSEXP);
    Rcpp::traits::input_parameter< int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma1(sigma1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::umat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double& >::type alpha_B(alpha_BSEXP);
    Rcpp::traits::input_parameter< double& >::type beta_B(beta_BSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type alpha_0(alpha_0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type alpha_1(alpha_1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_1(beta_1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_sampler(Y, Sd_hat, N, B, sigma, sigma1, sigma0, p, A, Z, alpha_B, beta_B, alpha_0, alpha_1, beta_0, beta_1, a, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FLOWMR_gibbs_sampler_with_corr", (DL_FUNC) &_FLOWMR_gibbs_sampler_with_corr, 21},
    {"_FLOWMR_gibbs_sampler", (DL_FUNC) &_FLOWMR_gibbs_sampler, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_FLOWMR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
