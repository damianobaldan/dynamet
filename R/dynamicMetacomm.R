#' Simulate Metacommunity Dynamics Over Multiple Epochs
#'
#' @param nEpochs The number of iterations/epochs to run.
#' @param ... All arguments passed to masterEqMetacomm. If an argument is a list,
#'            it must have length 1 (constant) or length nEpochs (dynamic).
#' @param init.comm Initial community state (optional).
#' @param updater A function that takes (current_comm, current_args) and
#'                returns an updated list of args for the next epoch.
#'
#' @returns A list of length nEpochs containing the resulting community matrix for each epoch.
#'
#' @importFrom utils modifyList
#'
#' @export
dynamicMetacomm <- function(nEpochs, ..., init.comm = NULL, updater = NULL) {

  #----------- 1. Input Validation ----------

  # Get default settings for the simulation
  defaults <- getSimulationDefaults()

  # Update with function input if something different from defaults was set
  args_list <- modifyList(defaults, list(...))

  # Remove init.comm from dots if user accidentally put it there
  args_list$init.comm <- NULL

  # Validate inputs
  validateDynamicInputs(
    args_list = args_list,
    nEpochs = nEpochs,
    init.comm = init.comm,
    updater = updater)

  # Initialize output list
  trajectory <- vector("list", nEpochs)

  # Initialize current community state
  current_comm <- init.comm

  #----------- 2. Iterate trough epochs ----------

  # Iterate trough epochs
  for (i in 1:nEpochs) {

    # Print current epoch
    if(args_list$verbose){
      message(sprintf("Running Epoch %d of %d...", i, nEpochs))
    }

    #----------- 2.1 Construct arguments for this specific epoch ----------

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
          return(arg)
        }
      })

    #----------- 2.2 Update init.comm and run the simulation ----------

    # Pass the current state
    epoch_args$init.comm <- current_comm

    # Run the core simulation
    current_comm <- do.call(masterEqMetacomm, epoch_args)

    # Save the result
    trajectory[[i]] <- current_comm

    #----------- 2.3 update parameters for the next simulation ----------

    # Update parameters if an updater function is provided
    if (!is.null(updater)) {

      # Calculate the new parameters based on the current community state and the current argument
      updated_params <- updater(current_comm, args_list)

      # Update the master list
      args_list <- utils::modifyList(args_list, updated_params)
    }
  }

  return(trajectory)
}



#' Validate inputs for dynamicMetacomm
#'
#' @param args_list A list of all arguments passed to dynamicMetacomm.
#' @param nEpochs The intended number of simulation epochs.
#'
#' @keywords internal
validateDynamicInputs <- function(args_list, nEpochs, init.comm, updater) {

  # Check if arguments that do not have a default are missing
  required <- c("Meta.pool", "Js", "M.migra", "m.pool")
  for (req in required) {
    if (is.null(args_list[[req]])) {
      stop(sprintf("Missing required argument: %s", req))
    }
  }

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

  # Validate updater function
  if (!is.null(updater)) {

    # Check if it is a function
    if (!is.function(updater)) {
      stop("'updater' must be a function.")
    }

    # Check the number of arguments
    # formals() returns a pairlist of the function's arguments
    updater_args <- names(formals(updater))

    if (length(updater_args) != 2) {
      stop(sprintf(
        "Invalid 'updater' function: It must accept exactly 2 arguments (current_comm, current_args), but it has %d.",
        length(updater_args)
      ))
    }
  }

}
