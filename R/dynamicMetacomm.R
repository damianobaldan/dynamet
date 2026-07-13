#' Simulate Metacommunity Dynamics Over Multiple Epochs
#'
#' @param nEpochs The number of iterations/epochs to run.
#' @param ... All arguments passed to masterEqMetacomm. If an argument is a list,
#'            it must have length 1 (constant) or length nEpochs (dynamic).
#' @param init.comm Initial community state (optional).
#'
#' @returns A list of length nEpochs containing the resulting community matrix for each epoch.
#'
#' @export
dynamicMetacomm <- function(nEpochs, ..., init.comm = NULL) {

  # Capture all arguments passed to the function
  args_list <- list(...)

  # Remove init.comm from dots if user accidentally put it there
  args_list$init.comm <- NULL

  # Validate
  validateDynamicInputs(args_list, nEpochs, init.comm)

  # Initialize output list
  trajectory <- vector("list", nEpochs)

  # Initialize current community state
  current_comm <- init.comm

  # Iterate trough epochs
  for (i in 1:nEpochs) {

    # Print current epoch
    message(sprintf("Running Epoch %d of %d...", i, nEpochs))

    # Construct arguments for this specific epoch:
    # apply function that checks if the argument is a list and
    # selects the appropriate index, otherwise returns the constant value.
    epoch_args <- lapply(
      args_list,
      function(arg) {
        # if argument is a list of length nEpochs take the i-th element, otherwise take the first element
        if (is.list(arg)) {
          index <- ifelse( length(arg) == nEpochs, i, 1)
          return(arg[[index]])
        # Otherwise simply return the argument
        } else {
          # Treat non-list as constant
          return(arg)
        }
      })

    # Pass the current state
    epoch_args$init.comm <- current_comm

    # Run the core simulation
    current_comm <- do.call(masterEqMetacomm, epoch_args)

    # Save the result
    trajectory[[i]] <- current_comm
  }

  return(trajectory)
}



#' Validate inputs for dynamicMetacomm
#'
#' @param args_list A list of all arguments passed to dynamicMetacomm.
#' @param nEpochs The intended number of simulation epochs.
#'
#' @keywords internal
validateDynamicInputs <- function(args_list, nEpochs, init.comm) {

  # Iterate trough the input arguments
  for (arg_name in names(args_list)) {

    # Get the current argument to be checked
    val <- args_list[[arg_name]]

    # If the user passed a list, we enforce that it is either
    # length 1 (constant) or length nEpochs (time-varying).
    if (is.list(val)) {
      if (length(val) != nEpochs && length(val) != 1) {
        stop(sprintf(
          "Validation Error: Argument '%s' is a list of length %d. It must be length 1 or %d.",
          arg_name, length(val), nEpochs
        ))
      }
    }
  }

  # Strictly validate init.comm
  if (!is.null(init.comm) && !is.matrix(init.comm)) {
    stop("'init.comm' must be NULL or a matrix representing the initial community state.")
  }


}
