#' Run Stochastic Replicates of dynamicMetacomm (Sequential or Parallel)
#'
#' @param nReplicates Integer. The number of independent simulation runs.
#' @param nEpochs Integer. Passed to dynamicMetacomm.
#' @param ... Arguments passed to dynamicMetacomm / masterEqMetacomm.
#' @param init.comm Matrix or NULL. Passed to dynamicMetacomm.
#' @param updater Function. Passed to dynamicMetacomm.
#' @param n_cores Integer. Number of CPU cores to use. If 1, runs sequentially without setting up a parallel backend.
#'
#' @returns A list of length nReplicates containing simulation trajectories.
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

  # Helper function to execute a single replicate call
  run_single_rep <- function(i) {
    do.call(dynamicMetacomm, c(
      list(
        nEpochs = nEpochs,
        init.comm = init.comm,
        updater = updater
      ),
      user_dots
    ))
  }

  # Sequential execution (n_cores == 1)
  if (n_cores == 1) {
    results <- lapply(seq_len(nReplicates), run_single_rep)

  } else {
    # Parallel execution (n_cores > 1)
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential), add = TRUE) # Guarantees cleanup on exit or error

    results <- future.apply::future_lapply(
      seq_len(nReplicates),
      run_single_rep,
      future.seed = TRUE
    )
  }

  return(results)
}
