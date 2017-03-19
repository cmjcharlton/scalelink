.onAttach <- function(...) {
  packageStartupMessage(paste(c("Scalelink: Record for 'probabilistic' linkage using a scaling procedure", capture.output(print(citation(package="Scalelink"))), "For examples on how to use this package type example(ScaleLink)"), collapse="\n"))
} 
