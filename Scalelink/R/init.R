.onAttach <- function(...) {
  packageStartupMessage(paste(c("Scalelink: Record for 'probabilistic' linkage using a scaling procedure", format(citation("Scalelink")), "\n", "For examples on how to use this package type example(Scalelink)"), collapse="\n"))
} 
