.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste("\nThis is LOBSTAHS version", packageVersion("LOBSTAHS"), "\n"))
}