#' Run Parallel Stochastic Replicates of dynamicMetacomm
#'
#' @description This function executes multiple independent simulations of
#' \code{dynamicMetacomm} in parallel using the future framework.
#'
#' @param nReplicates Integer. The number of independent simulation runs to perform.
#' @param nEpochs Integer. Passed to dynamicMetacomm.
#' @param ... Arguments passed to masterEqMetacomm (via dynamicMetacomm).
#' @param init.comm Matrix or NULL. Passed to dynamicMetacomm.
#' @param n_cores Integer. Number of CPU cores to use. Defaults to available cores - 1.
#'
#' @returns A list of length nReplicates, where each element is the resulting
#'          trajectory list from one simulation run.
#'
#' @import future
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
#'
#' @export
replicateDynamicMetacomm <- function(nReplicates, nEpochs, ..., init.comm = NULL, n_cores = parallel::detectCores() - 1) {

  # 1. Setup parallel backend
  # We use the explicitly namespaced calls here
  future::plan(future::multisession, workers = n_cores)

  # 2. Capture all ... arguments to pass them through
  sim_args <- list(...)

  # 3. Execute replicates
  # future_lapply is imported from future.apply
  results <- future.apply::future_lapply(seq_len(nReplicates), function(i) {

    # We call dynamicMetacomm with the captured parameters
    do.call(dynamicMetacomm, c(
      list(nEpochs = nEpochs, init.comm = init.comm),
      sim_args
    ))

  }, future.seed = TRUE) # Ensures stochastic independence across cores

  # 4. Clean up backend to return to normal
  future::plan(future::sequential)

  return(results)
}
