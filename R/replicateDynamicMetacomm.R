#' Run Stochastic Replicates of dynamicMetacomm (Sequential or Parallel)
#'
#' Executes multiple independent simulation replicates of \code{dynamicMetacomm}.
#' Supports both unipartite (standard) and bipartite mutualistic metacommunity
#' models in sequential or parallel processing modes.
#'
#' @param nReplicates Integer. The number of independent simulation runs.
#' @param nEpochs Integer. Passed to \code{dynamicMetacomm}.
#' @param ... Arguments passed to \code{dynamicMetacomm}, \code{masterEqMetacomm},
#'   or \code{masterEqMetacommMutualistic}.
#' @param init.comm Starting community state. For standard dynamics, a numeric matrix
#'   or NULL. For mutualistic dynamics, a list containing \code{Meta.animal} and
#'   \code{Meta.plant} matrices. Alternatively, a list of length \code{nReplicates}
#'   containing distinct initial conditions for each replicate.
#' @param updater Function. Optional parameter updater function passed to \code{dynamicMetacomm}.
#' @param n_cores Integer. Number of CPU cores to use. If 1, runs sequentially
#'   without setting up a parallel backend.
#'
#' @return A list of length \code{nReplicates}, where each element contains the
#'   epoch trajectory output from \code{dynamicMetacomm}.
#'
#' @import future
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
#' @export
replicateDynamicMetacomm <- function(nReplicates,
                                     nEpochs,
                                     ...,
                                     init.comm = NULL,
                                     updater = NULL,
                                     n_cores = max(1, parallel::detectCores() - 1)) {

  # Capture user dots cleanly
  user_dots <- list(...)

  # Check if init.comm is a single mutualistic initial state: list(Meta.animal=..., Meta.plant=...)
  is_mut_single_init <- is.list(init.comm) && all(c("Meta.animal", "Meta.plant") %in% names(init.comm))

  # Check if init.comm is a list containing separate initial conditions per replicate
  is_rep_specific_init <- is.list(init.comm) && length(init.comm) == nReplicates && !is_mut_single_init

  # Helper function to execute a single replicate call
  run_single_rep <- function(i) {
    # Extract replicate-specific initial state if provided
    rep_init <- if (is_rep_specific_init) init.comm[[i]] else init.comm

    do.call(dynamicMetacomm, c(
      list(
        nEpochs   = nEpochs,
        init.comm = rep_init,
        updater   = updater
      ),
      user_dots
    ))
  }

  # Sequential execution (n_cores == 1)
  if (n_cores == 1) {
    results <- lapply(seq_len(nReplicates), run_single_rep)

  # Parallel execution (n_cores > 1)
  } else {
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential), add = TRUE) # Guarantees backend cleanup on exit or error

    results <- future.apply::future_lapply(
      seq_len(nReplicates),
      run_single_rep,
      future.seed = TRUE
    )
  }

  return(results)
}
