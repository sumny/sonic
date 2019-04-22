register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1L)
  stopifnot(is.character(generic), length(generic) == 1L)
  stopifnot(is.character(class), length(class) == 1L)
  if(is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  } else {
    stopifnot(is.function(fun))
  }
  if(pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}

.onLoad <- function(libname, pkgname) {
  if(getRversion() < "3.6.0") {
    register_s3_method("sandwich", "estfun", "rm")
    register_s3_method("sandwich", "estfun", "twopl")
    register_s3_method("sandwich", "estfun", "threepl")
    register_s3_method("sandwich", "estfun", "threeplu")
    register_s3_method("sandwich", "estfun", "fourpl")
  }
  invisible()
}
