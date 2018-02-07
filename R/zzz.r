.onAttach <- function(lib, pkg) {
	packageStartupMessage("Welcome to GWmodel version 2.0-5.\n Note: This verision has been re-built with RcppArmadillo to improve its performance.\nThe new version of GWmodel 2.0-5 now is ready, and new functions have been incorporated:\n 1) gwr.psdm, GWR with parameter-specific distance metrics\n 2) gwr.bootstrap, Bootstrap GWR\n", appendLF = FALSE)
  }