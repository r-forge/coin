
.onLoad <- function(lib, pkg) {
    .Call("coin_init", PACKAGE = "coin")
    return(TRUE)
}

.onUnload <- function(libpath)
    library.dynam.unload("coin", libpath)
