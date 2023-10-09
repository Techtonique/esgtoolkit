


########## Interpolation

# Yield curve interpolation
ycInterpolation <- function(yM = NULL,
                            p = NULL,
                            matsin,
                            matsout,
                            method = c("NS", "SV", "SW", "HCSPL"),
                            typeres = c("rates", "prices"))
{
  if (missing(yM) &&
      missing(p))
    stop("market yield to maturities or prices must be provided")
  
  if (missing(matsin) ||
      missing(matsout))
    stop("Input and output maturities must be provided")
  
  if (max(matsout) > max(matsin))
    warning("Try extrapolation function instead. Check the output maturities.")
  
  method <- match.arg(method)
  
  typeres <- match.arg(typeres)
  
  if (is.null(yM) || missing(yM))
  {
    input_curve <- p
    # tZC_SW(UFR = UFR, T_UFR = T_UFR, u = matsin, p = input_curve, t = matsout, typeres = typeres)
    switch(
      method,
      NS = tZC_NS(
        p = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres
      ),
      SV = tZC_SV(
        p = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres
      ),
      SW = tZC_SW(
        UFR = 0.042,
        u = matsin,
        p = input_curve,
        t = matsout,
        typeres = typeres
      ),
      HCSPL = hermitecubicspline(
        p = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres
      )
    )
  }
  else {
    input_curve <- yM
    switch(
      method,
      NS = tZC_NS(
        yM = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres
      ),
      SV = tZC_SV(
        yM = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres
      ),
      SW = tZC_SW(
        UFR = 0.042,
        u = matsin,
        yM = input_curve,
        t = matsout,
        typeres = typeres
      ),
      HCSPL = hermitecubicspline(
        yM = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres
      )
    )
  }
}
ycInterpolation <- compiler::cmpfun(ycInterpolation)

########## Extrapolation

# Yield curve extrapolation wrapper
ycExtrapolation <- function(yM = NULL,
                            p = NULL,
                            matsin,
                            matsout,
                            method = c("NS", "SV", "SW"),
                            typeres = c("rates", "prices"),
                            UFR,
                            T_UFR = NULL)
{
  if (missing(yM) &&
      missing(p))
    stop("market yield to maturities or prices must be provided")
  if (missing(matsin) ||
      missing(matsout))
    stop("Input and output maturities must be provided")
  if (max(matsout) <= max(matsin))
    warning("max(matsout) <= max(matsin) : try the interpolation function instead")
  method <- match.arg(method)
  typeres <- match.arg(typeres)
  
  if (is.null(yM) || missing(yM))
  {
    input_curve <- p
    switch(
      method,
      NS = tZC_NS_extra(
        p = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres,
        UFR = UFR
      ),
      SV = tZC_SV_extra(
        p = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres,
        UFR = UFR
      ),
      SW = tZC_SW(
        UFR = UFR,
        T_UFR = T_UFR,
        u = matsin,
        p = input_curve,
        t = matsout,
        typeres = typeres
      )
    )
  }
  else {
    input_curve <- yM
    switch(
      method,
      NS = tZC_NS_extra(
        yM = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres,
        UFR = UFR
      ),
      SV = tZC_SV_extra(
        yM = input_curve,
        matsin = matsin,
        matsout = matsout,
        typeres = typeres,
        UFR = UFR
      ),
      SW = tZC_SW(
        UFR = UFR,
        T_UFR = T_UFR,
        u = matsin,
        yM = input_curve,
        t = matsout,
        typeres = typeres
      )
    )
  }
  
}
ycExtrapolation <- compiler::cmpfun(ycExtrapolation)