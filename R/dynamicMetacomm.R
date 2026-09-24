#' Simulate Metacommunity Dynamics Over Multiple Epochs
#'
#' @param nEpochs The number of iterations/epochs to run.
#' @param ... All arguments passed to \code{masterEqMetacomm} or
#'   \code{masterEqMetacommMutualistic}. If an argument is a list, it must have
#'   length 1 (constant) or length \code{nEpochs} (dynamic).
#' @param init.comm Initial community state (optional). For standard models, a matrix;
#'   for mutualistic models, a named list with \code{Meta.animal} and \code{Meta.plant} matrices.
#' @param updater A function that takes \code{(current_comm, current_args)} and
#'   returns an updated list of args for the next epoch.
#'
#' @return A list of length \code{nEpochs} containing the resulting community matrix
#'   (or list of matrices for mutualistic dynamics) for each epoch.
#'
#' @importFrom utils modifyList
#'
#' @export
dynamicMetacomm <- function(nEpochs, ..., init.comm = NULL, updater = NULL) {

  #----------- 1. Detect Simulation Mode & Setup Inputs ----------

  dots <- list(...)

  # Identify if mutualistic arguments are present
  mut_indicators <- c("Meta.pool.animal", "Meta.pool.plant", "Js.animal", "Js.plant",
                      "M.migra.animal", "M.migra.plant", "m.pool.animal", "m.pool.plant",
                      "alpha.animal.plant")
  is_mutualistic <- any(mut_indicators %in% names(dots))

  # Get default simulation settings
  defaults <- getSimulationDefaults()

  # Combine defaults with provided dot arguments
  args_list <- modifyList(defaults, dots)

  # Remove init.comm variants from args_list to prevent list indexing conflicts
  args_list$init.comm <- NULL
  args_list$init.comm.animal <- NULL
  args_list$init.comm.plant <- NULL

  # Initialize current community state
  if (is_mutualistic) {
    if (is.list(init.comm) && all(c("Meta.animal", "Meta.plant") %in% names(init.comm))) {
      current_comm <- init.comm
    } else if (!is.null(dots$init.comm.animal) || !is.null(dots$init.comm.plant)) {
      current_comm <- list(
        Meta.animal = dots$init.comm.animal,
        Meta.plant  = dots$init.comm.plant
      )
    } else {
      current_comm <- init.comm # NULL or user-provided list
    }
  } else {
    current_comm <- init.comm
  }

  # Validate inputs
  validateDynamicInputs(
    args_list      = args_list,
    nEpochs        = nEpochs,
    init.comm      = init.comm,
    updater        = updater,
    is_mutualistic = is_mutualistic,
    dots           = dots
  )

  # Initialize trajectory output list
  trajectory <- vector("list", nEpochs)

  #----------- 2. Iterate Through Epochs ----------

  for (i in 1:nEpochs) {

    if (isTRUE(args_list$verbose)) {
      message(sprintf("Running Epoch %d of %d (%s)...",
                      i, nEpochs, if (is_mutualistic) "Mutualistic" else "Standard"))
    }

    #----------- 2.1 Construct arguments for this specific epoch ----------
    epoch_args <- lapply(args_list, function(arg) {
      if (is.list(arg)) {
        index <- if (length(arg) == nEpochs) i else 1
        return(arg[[index]])
      } else {
        return(arg)
      }
    })

    #----------- 2.2 Attach init.comm and run simulation ----------
    if (is_mutualistic) {
      if (!is.null(current_comm)) {
        epoch_args$init.comm.animal <- current_comm$Meta.animal
        epoch_args$init.comm.plant  <- current_comm$Meta.plant
      }
      current_comm <- do.call(masterEqMutualistic, epoch_args)
    } else {
      if (!is.null(current_comm)) {
        epoch_args$init.comm <- current_comm
      }
      current_comm <- do.call(masterEqMetacomm, epoch_args)
    }

    # Save output for current epoch
    trajectory[[i]] <- current_comm

    #----------- 2.3 Parameter feedback updates (updater function) ----------
    if (!is.null(updater)) {
      updated_params <- updater(current_comm, epoch_args)
      args_list <- utils::modifyList(args_list, updated_params)
    }
  }

  return(trajectory)
}


#' Validate inputs for dynamicMetacomm
#'
#' @param args_list A list of all arguments passed to dynamicMetacomm.
#' @param nEpochs The intended number of simulation epochs.
#' @param init.comm Initial community input object.
#' @param updater Feedback updater function.
#' @param is_mutualistic Logical flag indicating mutualistic simulation mode.
#' @param dots Original dot arguments list.
#'
#' @keywords internal
validateDynamicInputs <- function(args_list, nEpochs, init.comm, updater, is_mutualistic, dots) {

  # Check required arguments based on simulation mode
  if (is_mutualistic) {
    required <- c("Meta.pool.animal", "Meta.pool.plant",
                  "Js.animal", "Js.plant",
                  "M.migra.animal", "M.migra.plant",
                  "m.pool.animal", "m.pool.plant")
  } else {
    required <- c("Meta.pool", "Js", "M.migra", "m.pool")
  }

  for (req in required) {
    if (is.null(args_list[[req]])) {
      stop(sprintf("Missing required argument for %s mode: %s",
                   if (is_mutualistic) "mutualistic" else "standard", req))
    }
  }

  # Validate length of time-varying list arguments
  for (arg_name in names(args_list)) {
    val <- args_list[[arg_name]]
    if (is.list(val)) {
      if (length(val) != nEpochs && length(val) != 1) {
        stop(sprintf(
          "Validation Error: Argument '%s' is a list of length %d. It must be length 1 or %d.",
          arg_name, length(val), nEpochs
        ))
      }
    }
  }

  # Validate initial community structure
  if (is_mutualistic) {
    if (!is.null(init.comm)) {
      if (is.list(init.comm)) {
        if (!all(c("Meta.animal", "Meta.plant") %in% names(init.comm)) ||
            !is.matrix(init.comm$Meta.animal) || !is.matrix(init.comm$Meta.plant)) {
          stop("'init.comm' for mutualistic dynamics must be NULL or a list with matrix elements 'Meta.animal' and 'Meta.plant'.")
        }
      } else if (!is.matrix(init.comm)) {
        stop("'init.comm' for mutualistic dynamics must be NULL or a list containing 'Meta.animal' and 'Meta.plant' matrices.")
      }
    }
    if (!is.null(dots$init.comm.animal) && !is.matrix(dots$init.comm.animal)) {
      stop("'init.comm.animal' must be a numeric matrix.")
    }
    if (!is.null(dots$init.comm.plant) && !is.matrix(dots$init.comm.plant)) {
      stop("'init.comm.plant' must be a numeric matrix.")
    }
  } else {
    if (!is.null(init.comm) && !is.matrix(init.comm)) {
      stop("'init.comm' must be NULL or a matrix representing the initial community state.")
    }
  }

  # Validate updater function signature
  if (!is.null(updater)) {
    if (!is.function(updater)) {
      stop("'updater' must be a function.")
    }

    updater_args <- names(formals(updater))
    if (length(updater_args) != 2) {
      stop(sprintf(
        "Invalid 'updater' function: It must accept exactly 2 arguments (current_comm, current_args), but it has %d.",
        length(updater_args)
      ))
    }
  }
}
