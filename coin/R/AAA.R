
if (FALSE) {
.First.lib <- function(libname, pkgname) {
    if (!require("survival")) 
        stop("cannot load ", sQuote("survival"))
    if (!require("mvtnorm")) 
        stop("cannot load ", sQuote("mvtnorm"))
    library.dynam("coin", pkgname, libname)
    .Call("coin_init", PACKAGE = "coin")
}

} else {

.onLoad <- function(lib, pkg) {
    if (!require("survival")) 
        stop("cannot load ", sQuote("survival"))
    if (!require("mvtnorm")) 
        stop("cannot load ", sQuote("mvtnorm"))
    .Call("coin_init", PACKAGE = "coin")
    return(TRUE)
}
}
