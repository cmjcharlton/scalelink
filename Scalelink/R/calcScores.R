#' Calculates linking scores for a file of interest and linkage data file.
#' 
#' This function calculates a score from two files, the file of interest (FOI)
#' and linkage data file (LDF).
#' 
#' @param FOI A \code{\link[base]{data.frame}} object, \code{\link[base]{matrix}} or 
#' \code{\link[base]{vector}} to be used as the file of interest. This must contain
#' only the variables of interest, and these must be in the same order as the LDF.
#'
#' @param LDF A \code{\link[base]{data.frame}} object, \code{\link[base]{matrix}} or 
#' \code{\link[base]{vector}} to be used as the linkage data file. This must contain
#' only the variables of interest, and these must be in the same order as the FOI.
#'
#' @param missing.value Value used to represent missing data; Defaults to NA
#'
#' @param min.parallelblocksize The minimum block size when splitting up the data accross
#' processors. You may wish to change this to optimise the allocation of processors.
#' see (\url{https://rcppcore.github.io/RcppParallel/#tuning}).
#'
#' @param output.varnames Labels to apply to function output; Defaults to column names
#' of the FOI \code{\link[base]{data.frame}}
#'
#' @param debug Boolean indicating whether to output additional debugging information
#' 
#' @return A list containing: An numeric \code{\link[base]{vector}} of scores, one for
#' each of the identifiers of interest and a \code{\link[base]{matrix}} containing A*.
#' 
#' @author Goldstein H., and Charlton, C.M.J., (2017) Centre for Multilevel Modelling, University of Bristol.
#' 
#' @export
calcScores <- function(FOI, LDF, missing.value=NA, min.parallelblocksize=1, output.varnames=NULL, debug=FALSE) {
  NFOI <- nrow(FOI)
  NLDF <- nrow(LDF)
  NIDENT <- ncol(FOI)

  if (is.null(output.varnames)) {
    if (!all(colnames(FOI) == colnames(LDF))) {
      stop("variable names in FOI and LDF do not match")
    }
    output.varnames <- colnames(FOI)
  } else {
    colnames(FOI) <- NULL
    colnames(LDF) <- NULL
  }

  if (NIDENT != ncol(LDF)) {
    stop("Number of identifiers in file of interest is different to number in linking file")
  }

  if (length(output.varnames) != ncol(LDF)) {
    stop("Incorrect number of variable names specified in output.varnames")
  }

  if (min.parallelblocksize < 0 || min.parallelblocksize > nrow(FOI)) {
    stop("Invalid value for min.parallelblocksize")
  }

  starttime <- proc.time()

  if (!is.na(missing.value)) {
    FOI[FOI==missing.value] <- NA
    LDF[LDF==missing.value] <- NA
  }
  
  # randomly insert for missing values - first LDF then FOI
  LDFCOMP <- complete.cases(LDF)
  LDFNEW <- LDF[LDFCOMP, , drop=FALSE]
  LDF[!LDFCOMP, ] <- LDFNEW[sample(nrow(LDFNEW), size=(nrow(LDF) - nrow(LDFNEW)), replace=TRUE), ]

  FOICOMP <- complete.cases(FOI)
  FOINEW <- FOI[FOICOMP, , drop=FALSE]
  FOI[!FOICOMP, ] <- FOINEW[sample(nrow(FOINEW), size=(nrow(FOI) - nrow(FOINEW)), replace=TRUE), ]
  
  #  Accumulates matrix values. Note: assumes all identifiers have just 2 categories.
  ASTAR <- buildAstar(data.frame(FOI), data.frame(LDF), min.parallelblocksize, debug)
  if (isTRUE(debug)) {
    print("A*:")
    print(ASTAR)
  }
  B <- matrix(data=c(rep(0, (NIDENT * 2) + 1), 1 / NIDENT) * NFOI * NLDF, ncol=1)
  if (isTRUE(debug)) {
    print("B:")
    print(B)
  }
  SCORES <- solve(ASTAR, B)

  if (isTRUE(debug)) {
    print("SCORES")
    print(SCORES)
  }

  ADD <- rep(SCORES[((1:NIDENT) * 2) - 1, 1], each=2)
  SCORES[1:(NIDENT*2), 1] <- (SCORES[1:(NIDENT*2), 1] - ADD) * 100

  if (isTRUE(debug)) {
    print("Adjusted SCORES")
    print(SCORES)
  }
 
  #print(paste("Number_of_primary_file_records", NFOI))

  OUTSCORES <- SCORES[1:NIDENT*2, ]
  # Scale scores
  OUTSCORES <- OUTSCORES / sum(OUTSCORES) * 100
  names(OUTSCORES) <- output.varnames

  endtime <- proc.time() - starttime
  cat(paste("Elapsed time:", endtime[3], "\n"))

  list(scores=OUTSCORES, astar=ASTAR)
}