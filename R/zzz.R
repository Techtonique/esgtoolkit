.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n \n This is version ",
    utils::packageVersion(pkgname),
    " of ",
    pkgname,
    ". Starting with 1.0.0, it will be renamed as: 'esgtoolkit' (to finally remove all my active packages from CRAN) \n \n "
  )
} 