# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rnormESGcpp <- function(N, M) {
    .Call('_esgtoolkit_rnormESGcpp', PACKAGE = 'esgtoolkit', N, M)
}

rOUESGcpp <- function(N, horizon, Delta, x0, theta, eps) {
    .Call('_esgtoolkit_rOUESGcpp', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta, eps)
}

rOUESGcppexact <- function(N, horizon, Delta, x0, theta, eps) {
    .Call('_esgtoolkit_rOUESGcppexact', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta, eps)
}

rCIRESGcpp <- function(N, horizon, Delta, x0, theta, eps) {
    .Call('_esgtoolkit_rCIRESGcpp', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta, eps)
}

rCIRESGcppexact <- function(N, horizon, Delta, x0, theta, eps) {
    .Call('_esgtoolkit_rCIRESGcppexact', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta, eps)
}

rGBMESGcpp <- function(N, horizon, Delta, x0, theta1, theta2, eps) {
    .Call('_esgtoolkit_rGBMESGcpp', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta1, theta2, eps)
}

rGBMjumpsnormESGcpp <- function(N, horizon, Delta, x0, theta1, theta2, lambda, mu, sigma, eps) {
    .Call('_esgtoolkit_rGBMjumpsnormESGcpp', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta1, theta2, lambda, mu, sigma, eps)
}

rGBMjumpskouESGcpp <- function(N, horizon, Delta, x0, theta1, theta2, lambda, eta_up, eta_down, p, eps) {
    .Call('_esgtoolkit_rGBMjumpskouESGcpp', PACKAGE = 'esgtoolkit', N, horizon, Delta, x0, theta1, theta2, lambda, eta_up, eta_down, p, eps)
}

TAGcorecpp <- function(sim, sj_down, sj_up, n, p) {
    .Call('_esgtoolkit_TAGcorecpp', PACKAGE = 'esgtoolkit', sim, sj_down, sj_up, n, p)
}

