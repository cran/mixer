.onLoad <- function(libname, pkgname) {
  #packageStartupMessage("\n")
  #packageStartupMessage("      'mixer' package\n"  )
  #packageStartupMessage("      SSB group (http://ssbgroup.fr)\n")
  #packageStartupMessage("      mixer page (http://ssbgroup.fr/mixnet/mixer.html)\n")
  #packageStartupMessage("\n")
    
  setSeed(  Sys.getpid() )
}


