
call_libcoin <- function(routine, ...)
    .Call(routine, ..., PACKAGE = "libcoin")
