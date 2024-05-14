.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n \n This is version ",
    utils::packageVersion(pkgname),
    " of ",
    pkgname,
    ". Starting with 1.0.0, package renamed as: 'esgtoolkit' (lowercase) \n \n "
  )
} 