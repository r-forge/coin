
.First.lib <- function(libname, pkgname) {
    if (!require("survival")) 
        stop("cannot load ", sQuote("survival"))
    if (!require("mvtnorm")) 
        stop("cannot load ", sQuote("mvtnorm"))
    library.dynam("coin", pkgname, libname)
    .Call("coin_init", PACKAGE = "coin")
}
