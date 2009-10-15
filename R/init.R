.onLoad <- function(libname, pkgname) {
  cat("\n")
  cat("      'mixer' package\n"  )
  cat("      SSB group (http://ssbgroup.fr)\n")
  cat("      mixer page (http://ssbgroup.fr/mixnet/mixer.html)\n")
  cat("\n")
    
  setSeed(  Sys.getpid() )
}


