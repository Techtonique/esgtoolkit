// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rnormESGcpp
NumericMatrix rnormESGcpp(const int N, const int M);
RcppExport SEXP _esgtoolkit_rnormESGcpp(SEXP NSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(rnormESGcpp(N, M));
    return rcpp_result_gen;
END_RCPP
}
// rOUESGcpp
NumericMatrix rOUESGcpp(const int N, const int horizon, const double Delta, const double x0, NumericVector theta, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rOUESGcpp(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP thetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rOUESGcpp(N, horizon, Delta, x0, theta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rOUESGcppexact
NumericMatrix rOUESGcppexact(const int N, const int horizon, const double Delta, const double x0, NumericVector theta, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rOUESGcppexact(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP thetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rOUESGcppexact(N, horizon, Delta, x0, theta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rCIRESGcpp
NumericMatrix rCIRESGcpp(const int N, const int horizon, const double Delta, const double x0, NumericVector theta, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rCIRESGcpp(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP thetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rCIRESGcpp(N, horizon, Delta, x0, theta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rCIRESGcppexact
NumericMatrix rCIRESGcppexact(const int N, const int horizon, const double Delta, const double x0, NumericVector theta, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rCIRESGcppexact(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP thetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rCIRESGcppexact(N, horizon, Delta, x0, theta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rGBMESGcpp
NumericMatrix rGBMESGcpp(const int N, const int horizon, const double Delta, const double x0, NumericMatrix theta1, NumericMatrix theta2, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rGBMESGcpp(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP theta1SEXP, SEXP theta2SEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rGBMESGcpp(N, horizon, Delta, x0, theta1, theta2, eps));
    return rcpp_result_gen;
END_RCPP
}
// rGBMjumpsnormESGcpp
NumericMatrix rGBMjumpsnormESGcpp(const int N, const int horizon, const double Delta, const double x0, NumericMatrix theta1, NumericMatrix theta2, const double lambda, const double mu, const double sigma, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rGBMjumpsnormESGcpp(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP theta1SEXP, SEXP theta2SEXP, SEXP lambdaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rGBMjumpsnormESGcpp(N, horizon, Delta, x0, theta1, theta2, lambda, mu, sigma, eps));
    return rcpp_result_gen;
END_RCPP
}
// rGBMjumpskouESGcpp
NumericMatrix rGBMjumpskouESGcpp(const int N, const int horizon, const double Delta, const double x0, NumericMatrix theta1, NumericMatrix theta2, const double lambda, const double eta_up, const double eta_down, const double p, NumericMatrix eps);
RcppExport SEXP _esgtoolkit_rGBMjumpskouESGcpp(SEXP NSEXP, SEXP horizonSEXP, SEXP DeltaSEXP, SEXP x0SEXP, SEXP theta1SEXP, SEXP theta2SEXP, SEXP lambdaSEXP, SEXP eta_upSEXP, SEXP eta_downSEXP, SEXP pSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const double >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type eta_up(eta_upSEXP);
    Rcpp::traits::input_parameter< const double >::type eta_down(eta_downSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rGBMjumpskouESGcpp(N, horizon, Delta, x0, theta1, theta2, lambda, eta_up, eta_down, p, eps));
    return rcpp_result_gen;
END_RCPP
}
// TAGcorecpp
NumericVector TAGcorecpp(NumericVector sim, NumericVector sj_down, NumericVector sj_up, const int n, const int p);
RcppExport SEXP _esgtoolkit_TAGcorecpp(SEXP simSEXP, SEXP sj_downSEXP, SEXP sj_upSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sim(simSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sj_down(sj_downSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sj_up(sj_upSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(TAGcorecpp(sim, sj_down, sj_up, n, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_esgtoolkit_rnormESGcpp", (DL_FUNC) &_esgtoolkit_rnormESGcpp, 2},
    {"_esgtoolkit_rOUESGcpp", (DL_FUNC) &_esgtoolkit_rOUESGcpp, 6},
    {"_esgtoolkit_rOUESGcppexact", (DL_FUNC) &_esgtoolkit_rOUESGcppexact, 6},
    {"_esgtoolkit_rCIRESGcpp", (DL_FUNC) &_esgtoolkit_rCIRESGcpp, 6},
    {"_esgtoolkit_rCIRESGcppexact", (DL_FUNC) &_esgtoolkit_rCIRESGcppexact, 6},
    {"_esgtoolkit_rGBMESGcpp", (DL_FUNC) &_esgtoolkit_rGBMESGcpp, 7},
    {"_esgtoolkit_rGBMjumpsnormESGcpp", (DL_FUNC) &_esgtoolkit_rGBMjumpsnormESGcpp, 10},
    {"_esgtoolkit_rGBMjumpskouESGcpp", (DL_FUNC) &_esgtoolkit_rGBMjumpskouESGcpp, 11},
    {"_esgtoolkit_TAGcorecpp", (DL_FUNC) &_esgtoolkit_TAGcorecpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_esgtoolkit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
