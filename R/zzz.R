.onLoad <- function(libname, pkgname) {
  # Register the parsnip engine if parsnip is already loaded, and otherwise
  # when it is loaded later. parsnip is only suggested, so it is not loaded
  # here.
  if (isNamespaceLoaded("parsnip")) {
    try(.fcr_register_count_reg(), silent = TRUE)
  } else {
    setHook(packageEvent("parsnip", "onLoad"), function(...) {
      try(.fcr_register_count_reg(), silent = TRUE)
    })
  }
  invisible()
}
