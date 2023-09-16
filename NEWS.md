# version 0.6.0

- include hermite cubic spline and Smith-Wilson methods for curve interpolation

# version 0.5.0

- Calculate returns or log-returns for multivariate time series with `calculatereturns`
- Specify `start` (the time of the first observation, see `stats::ts` for examples) for simulations in `simdiff` and `simshocks`

# version 0.4.0

- use spline interpolation (stats::spline) for forward rates 


# version 0.3.0

- Use [VineCopula](http://tnagler.github.io/VineCopula/) package instead of CDVine (archived) for dependency simulation (means there's also now VineCopula's R-Vine Copulas simulation)
- Remove roxygen2 comments
- Create website with pkgdown, including docs --> https://techtonique.github.io/ESGtoolkit/


# version 0.2.0

- refactor files 
- add reproductibility seeds
- update vignette
- update LICENSE

# version 0.1.0

- Initial version