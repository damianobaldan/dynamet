#' Simulate Metacommunity Dynamics via Coalescent Assembly and lottery Phases
#'
#' @description
#' This function models metacommunity assembly and demographic turnover across a
#' network of local communities. It initiates patches using a coalescent
#' process up to local carrying capacities (\code{Js}). Then, an optional neutral/niche
#' lottery dynamic is used to simulate metacommunities. The simulation incorporates spatial migration, regional pool
#' immigration, self-recruitment, species-specific dispersal constraints (\code{d.spp}), environmental filtering (\code{FF}),
#' interspecific competition (\code{alpha}), and temperature-dependent mortality
#' scaled via the Arrhenius equation.
#'
#' @param Meta.pool A numeric vector of length S representing the relative abundances
#'   or probabilities of species within the regional species pool.
#' @param Js A numeric vector of length C setting the local carrying capacity (total
#'   individual slots) for each community patch.
#' @param M.migra A square numeric matrix of dimensions C x C establishing spatial
#'   migration connectivity and dispersal probabilities between patches. Cannot be NULL.
#' @param m.pool A single numeric value between 0 and 1 defining the probability of
#'   recruitment originating from the global regional pool rather than local/neighboring sources.
#' @param d.spp An optional numeric vector of length S dictating species-specific dispersal traits.
#'   If \code{NULL} (default), it is treated as a vector of 1s (equal dispersal ability for all species).
#' @param FF An optional numeric matrix of dimensions S x C representing local environmental
#'   filtering filters. Individual coefficients must range between 0 and 1, where 1 indicates
#'   perfect environmental match (no filtering penalty). If \code{NULL} (default), no environmental filtering is applied.
#' @param alpha An optional square numeric matrix of dimensions S x S defining interspecific
#'   competition coefficients between all pairs of species. If \code{NULL} (default), no interspecific competition is applied.
#' @param init.comm An optional numeric matrix of dimensions S x C representing the custom
#'   starting abundance counts for all species across patches. Required if \code{coalescence = FALSE}.
#' @param id.fixed An optional numeric vector containing indices of communities whose
#'   compositions are static and locked during simulation steps.
#' @param comm.fixed An optional numeric vector of length S outlining the fixed relative
#'   abundance profile assigned to the patches listed in \code{id.fixed}. These relative
#'   abundances are automatically scaled to the specific local carrying capacity (\code{Js})
#'   of each fixed patch.
#' @param prop.dead.by.it A numeric fraction (0, 1) defining the baseline mortality
#'   turnover rate for the coldest patch.
#' @param Ea A numeric value representing the activation energy (in eV) utilized to
#'   scale temperature-driven mortality.
#' @param Ts A numeric vector of length C specifying the local patch ambient temperatures
#'   expressed strictly in Kelvin units (> 0).
#' @param m.temp An optional numeric vector of length C, or a matrix of dimensions S x C,
#'   dictating community memory weight parameters (e.g., local seed banks or persistent vegetative state).
#' @param lottery A single logical value. If \code{TRUE}, triggers the demographic turnover
#'   iterative lottery phase after initial community assembly.
#' @param nIterations A single numeric integer specifying the total timeline steps/iterations to
#'   run inside the lottery loop.
#' @param verbose A single logical value. If \code{TRUE} (default), outputs live processing
#'   milestones and loop updates to the console.
#'
#' @returns A numeric matrix of dimensions S x C where cells hold the absolute abundance
#'   counts of individuals for each species across all simulated community patches.
#'
#' @author Matias Arim & Ana Borthagaray
#'
#' @export
#'
#' @importFrom stats rmultinom
#'
masterEqMetacomm <- function(Meta.pool, Js, M.migra, m.pool, d.spp = NULL,
                             FF = NULL, alpha = NULL, init.comm = NULL,
                             id.fixed = NULL, comm.fixed = NULL,
                             prop.dead.by.it = 0.05, Ea = 1e-5, Ts = 293.15, m.temp = 0,
                             lottery = TRUE, nIterations = 100, verbose = TRUE) {


  #----------- 1. VALIDATE INPUTS -----------
  validateMetaInputs(
    Meta.pool       = Meta.pool,
    Js              = Js,
    M.migra         = M.migra,
    m.pool          = m.pool,
    d.spp           = d.spp,
    FF              = FF,
    alpha           = alpha,
    init.comm       = init.comm,
    id.fixed        = id.fixed,
    comm.fixed      = comm.fixed,
    prop.dead.by.it = prop.dead.by.it,
    Ea              = Ea,
    Ts              = Ts,
    m.temp          = m.temp,
    lottery         = lottery,
    nIterations     = nIterations
  )

  #----------- 2. INITIALIZATION AND DATA NORMALIZATION -----------

  # Structural Dimension References
  S <- length(Meta.pool)
  C <- length(Js)

  # Handle d.spp default value
  if (is.null(d.spp)) {
    d.spp <- rep(1, S)
  }

  # Normalize vectors to represent relative probabilities summing to 1
  d.spp      <- d.spp / sum(d.spp)
  Meta.pool  <- Meta.pool / sum(Meta.pool)

  # Normalize comm.fixed if it is not NULL
  if (!is.null(comm.fixed)) {
    if (is.matrix(comm.fixed)) {
      # If matrix: sweep through and normalize each column independently
      comm.fixed <- sweep(comm.fixed, 2, colSums(comm.fixed), FUN = "/")
    } else {
      # If vector: normalize normally
      comm.fixed <- comm.fixed / sum(comm.fixed)
    }
  }

  # If FF is NULL, initialize it with 1s (No filtering effect)
  if (is.null(FF)) {
    FF <- matrix(1, nrow = S, ncol = C)
  }

  # Initialize community memory tracking (m.temp) as a strict S x C matrix
  if (is.null(m.temp)) {
    m.temp <- matrix(0, nrow = S, ncol = C)
  } else if (!is.matrix(m.temp)) {
    if (length(m.temp) == 1) {
      m.temp <- matrix(m.temp, nrow = S, ncol = C)
    } else if (length(m.temp) == C) {
      m.temp <- matrix(rep(m.temp, each = S), nrow = S, ncol = C)
    }
  }


  #----------- 3. THERMAL DEPENDENCE CALCULATIONS -----------

  # Estimate standard constant based on the minimum temperature patch to anchor baseline mortality
  min.dead.Tmin <- prop.dead.by.it / exp(-Ea / (min(Ts) * 8.62e-5))

  # Scale individual community mortality rates based on local temperature
  prop.dead.by.comm <- min.dead.Tmin * exp(-Ea / (Ts * 8.62e-5))

  # Prevent simulation failure if temperature-driven mortality is unrealistically high
  if (max(prop.dead.by.comm) > 0.95) {
    return("Error: Higher temperature replaces >100% of individuals. prop.dead.by.it is likely too large.")
  }

  # Determine integer count of deaths per community per iteration
  dead.by.it <- round(prop.dead.by.comm * Js, 0)
  dead.by.it <- ifelse(dead.by.it < 2, 2, dead.by.it) # Enforce a minimum floor of 2 deaths

  # Avoid fixed communities from undergoing standard mortality
  if (!is.null(id.fixed)) {
    dead.by.it[id.fixed] <- 0
  }
  max.dead.by.it <- max(dead.by.it)


  #----------- 4. COALESCENT ASSEMBLY PHASE -----------

  # If no initial community is provided, seed each community with exactly 1 individual based on regional pool and filters
  if (is.null(init.comm)) {
    Meta <- matrix(0, nrow = S, ncol = C)
    for (i in 1:ncol(M.migra)) {
      Meta[,i] <- stats::rmultinom(1, 1, Meta.pool * d.spp * FF[,i])
    }
  }

  # If initial community is provided, this becames the new metacommunity matrix
  if (!is.null(init.comm)) {
    Meta <- init.comm
  }

  # Loop over communities to perform the coalescent assembly
  for (ii in 1:max(Js)) {

    # Target communities that are under capacity AND whose current individual count is below ii
    id.j <- which(colSums(Meta) < Js & colSums(Meta) < ii)

    # Exclude fixed communities from receiving random individuals (they are managed below)
    if (!is.null(id.fixed)) {
      id.j <- setdiff(id.j, id.fixed)
    }

    # Check: if no communities need an individual in this step, skip to next iteration
    if (length(id.j) == 0){next}

    # Keep fixed communities scaled to current global abundance level
    if (!is.null(id.fixed)) {
      for (idx in seq_along(id.fixed)) {
        f_id <- id.fixed[idx]
        current_profile <- if (is.matrix(comm.fixed)) comm.fixed[, idx] else comm.fixed
        Meta[, f_id] <- current_profile * min(ii - 1, Js[f_id])
      }
    }

    # Calculate potential recruits based on local abundance and spatial migration
    Pool.neighbor <- (Meta %*% M.migra) * d.spp * FF

    # Calculate interspecific competition overlap (this is skipped if alpha is NULL)
    if (!is.null(alpha)) {
      overlap <- alpha %*% Meta
      col_sums_overlap <- colSums(overlap)
      col_sums_overlap[col_sums_overlap == 0] <- 1
      overlap <- sweep(overlap, 2, col_sums_overlap, FUN = "/")
      Pool.neighbor <- Pool.neighbor * (1 - overlap)
    }

    # If Pool.neighbor drops to 0, use regional pool layout as fallback
    col_sums_coalescent <- colSums(Pool.neighbor)
    if (any(col_sums_coalescent == 0)) {
      zero_cols <- which(col_sums_coalescent == 0)
      for (zc in zero_cols) {
        Pool.neighbor[, zc] <- Meta.pool
      }
    }

    # Assign new individuals to available spaces
    if (length(id.j) > 1) {
      new <- apply(Pool.neighbor[, id.j, drop = FALSE], 2, born, dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool)
      Meta[, id.j] <- Meta[, id.j] + new
    } else {
      Meta[, id.j] <- Meta[, id.j] + born(n = Pool.neighbor[, id.j], dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool)
    }

    if(verbose){ cat("coalescent construction in J:", ii, "of", max(Js), "\n") }
  }

  #----------- 5. LOTTERY DYNAMICS -----------
  if (lottery) {

    # Mirror state for community memory tracking
    Meta.lag <- Meta

    # Set fixed communities safely scaled to their local individual capacities (Js)
    if (!is.null(id.fixed)) {
      for (idx in seq_along(id.fixed)) {
        f_id <- id.fixed[idx]
        current_profile <- if (is.matrix(comm.fixed)) comm.fixed[, idx] else comm.fixed
        Meta[, f_id] <- round(current_profile * Js[f_id], 0)
      }
    }

    # Generate timeline checkpoints
    generations <- seq(1, nIterations, 1 / prop.dead.by.it)

    # Loop over iterations of the lottery dynamic
    for (iteration in 1:nIterations) {

      # Update community memory cache if at checkpoint
      if (iteration %in% generations) { Meta.lag <- Meta }

      # --- Death Sub-phase ---
      for (dead in 1:max.dead.by.it) {

        # Target patches experiencing active mortality
        id.dead <- which(dead.by.it >= dead)

        # Select and remove individuals
        if (length(id.dead) > 1) {
          Meta[, id.dead] <- Meta[, id.dead] - apply(Meta[, id.dead] * (1.001 - FF[, id.dead]), 2, FUN = change, change = 1)
        }
        if (length(id.dead) == 1) {
          Meta[, id.dead] <- Meta[, id.dead] - change(n = Meta[, id.dead] * (1.001 - FF[, id.dead]), change = 1)
        }
      }

      # --- Recruitment Sub-phase ---

      # Calculate potential recruits based on local abundance and spatial migration
      Pool.neighbor <- (Meta %*% M.migra) * d.spp * FF

      # Factor in community memory matrix element-wise
      Pool.neighbor <- Pool.neighbor * (1 - m.temp) + Meta.lag * (m.temp)

      # Re-calculate interspecific competition overlap (this is skipped if alpha is NULL)
      if (!is.null(alpha)) {
        overlap <- alpha %*% Meta
        col_sums_overlap <- colSums(overlap)
        col_sums_overlap[col_sums_overlap == 0] <- 1
        overlap <- sweep(overlap, 2, col_sums_overlap, FUN = "/")
        Pool.neighbor <- Pool.neighbor * (1 - overlap)
      }

      # Normalize Pool.neighbor safely by column
      col_sums <- colSums(Pool.neighbor)
      col_sums[col_sums == 0] <- 1
      Pool.norm <- sweep(Pool.neighbor, 2, col_sums, FUN = "/")

      Prob.mat  <- (1 - m.pool) * Pool.norm + m.pool * Meta.pool

      # Fill vacancies using multinomial sampling across local, neighboring, and regional pools
      for (i in seq_along(dead.by.it)) {
        Meta[, i] <- Meta[, i] + stats::rmultinom(1, size = dead.by.it[i], prob = Prob.mat[, i])
      }

      if(verbose){ cat("lottery iteration", iteration, "of", nIterations, "\n")}
    }
  }

  return(Meta)
}


#' Validate the inputs for masterEqMetacomm
#'
#' @inheritParams masterEqMetacomm
#' @keywords internal
#'
validateMetaInputs <- function(
    Meta.pool,
    Js,
    M.migra,
    m.pool,
    d.spp,
    FF = NULL,
    alpha = NULL,
    init.comm = NULL,
    id.fixed = NULL,
    comm.fixed = NULL,
    prop.dead.by.it,
    Ea,
    Ts,
    m.temp = NULL,
    lottery,
    nIterations) {

  # ----------------------------------------------------------------------------
  # 1. CLASS AND TYPE VALIDATIONS
  # ----------------------------------------------------------------------------
  if (!is.numeric(m.pool) || length(m.pool) != 1) stop("'m.pool' must be a single numeric value.")
  if (!is.logical(lottery) || length(lottery) != 1) stop("'lottery' must be a single logical value.")
  if (!is.numeric(nIterations) || length(nIterations) != 1) stop("'nIterations' must be a single numeric integer.")
  if (!is.numeric(prop.dead.by.it) || length(prop.dead.by.it) != 1) stop("'prop.dead.by.it' must be a single numeric fraction.")
  if (!is.numeric(Ea) || length(Ea) != 1) stop("'Ea' must be a single numeric value.")

  if (!is.numeric(Meta.pool)) stop("'Meta.pool' must be a numeric vector.")
  if (!is.numeric(d.spp))      stop("'d.spp' must be a numeric vector.")
  if (!is.numeric(Js))         stop("'Js' must be a numeric vector.")
  if (!is.numeric(Ts))         stop("'Ts' must be a numeric vector.")

  if (!is.null(FF) && (!is.matrix(FF) || !is.numeric(FF)))     stop("'FF' must be a numeric matrix.")
  if (!is.null(alpha) && (!is.matrix(alpha) || !is.numeric(alpha))) stop("'alpha' must be a numeric matrix.")

  # Configuration cross-dependence safety
  if (!is.null(id.fixed) && is.null(comm.fixed)) {
    stop("Configuration error: 'id.fixed' was provided, but 'comm.fixed' is NULL.")
  }
  if (is.null(id.fixed) && !is.null(comm.fixed)) {
    stop("Configuration error: 'comm.fixed' was provided, but 'id.fixed' is NULL.")
  }

  # Allow comm.fixed to be either a vector or a matrix
  if (!is.null(comm.fixed) && !is.numeric(comm.fixed)) {
    stop("'comm.fixed' must be a numeric vector or matrix.")
  }

  if (!is.null(M.migra) && (!is.matrix(M.migra) || !is.numeric(M.migra))) stop("'M.migra' must be a numeric matrix.")
  if (!is.null(id.fixed) && !is.numeric(id.fixed)) stop("'id.fixed' must be a numeric vector.")
  if (!is.null(m.temp) && !is.numeric(m.temp))     stop("'m.temp' must be a numeric vector or matrix.")

  # ----------------------------------------------------------------------------
  # 2. DIMENSION CONSISTENCY VALIDATIONS
  # ----------------------------------------------------------------------------
  S <- length(Meta.pool)
  C <- length(Js)

  # Species-dimension alignments (Rows)
  if (length(d.spp) != S) stop(sprintf("Dimension mismatch: 'd.spp' length (%d) must match 'Meta.pool' (%d).", length(d.spp), S))
  if (!is.null(FF) && nrow(FF) != S)       stop(sprintf("Dimension mismatch: Environmental filter matrix 'FF' must have %d rows (Species).", S))
  if (!is.null(alpha) && (nrow(alpha) != S || ncol(alpha) != S)) stop(sprintf("Dimension mismatch: Interspecific competition matrix 'alpha' must be a square matrix of %d x %d.", S, S))

  # Matrix/Vector dimension handler for comm.fixed
  if (!is.null(comm.fixed)) {
    if (is.matrix(comm.fixed)) {
      if (nrow(comm.fixed) != S || ncol(comm.fixed) != length(id.fixed)) {
        stop(sprintf("Dimension mismatch: 'comm.fixed' matrix must have dimensions S (%d) x length(id.fixed) (%d).", S, length(id.fixed)))
      }
    } else {
      if (length(comm.fixed) != S) {
        stop(sprintf("Dimension mismatch: 'comm.fixed' vector must be of length S (%d).", S))
      }
    }
  }

  # Community-dimension alignments (Columns)
  if (!is.null(FF) && ncol(FF) != C) stop(sprintf("Dimension mismatch: Environmental filter matrix 'FF' must have %d columns (Communities).", C))

  if (!is.null(M.migra)) {
    if (nrow(M.migra) != C || ncol(M.migra) != C) stop(sprintf("Dimension mismatch: Migration matrix 'M.migra' must be a square matrix of %d x %d.", C, C))
  } else {
    stop("Missing required input: 'M.migra' matrix cannot be NULL.")
  }

  # Handle m.temp dimension check gracefully if it is passed as a matrix
  if (!is.null(m.temp)) {
    if (is.matrix(m.temp)) {
      if (nrow(m.temp) != S || ncol(m.temp) != C) {
        stop(sprintf("Dimension mismatch: Matrix 'm.temp' must match community landscape dimensions (%d rows x %d columns).", S, C))
      }
    } else if (length(m.temp) > 1 && length(m.temp) != C) {
      stop(sprintf("Dimension mismatch: Vector 'm.temp' length (%d) must match the number of communities (%d).", length(m.temp), C))
    }
  }

  if (!is.null(id.fixed) && (any(id.fixed < 1) || any(id.fixed > C))) stop("'id.fixed' contains out-of-bounds community indices.")

  # ----------------------------------------------------------------------------
  # 3. VALUE AND BOUNDARY CONSTRAINTS (Includes Zero-Sum / NaN Protections)
  # ----------------------------------------------------------------------------
  # Missing data sweep
  all_inputs <- list(Meta.pool, d.spp, m.pool, Js, M.migra, nIterations, prop.dead.by.it, Ea, Ts)
  if (any(sapply(all_inputs, function(x) any(is.na(x))))) stop("Value error: Missing values (NA/NaN) detected in primary inputs.")
  if (!is.null(FF) && any(is.na(FF))) stop("Value error: Missing values (NA/NaN) detected in matrix 'FF'.")
  if (!is.null(alpha) && any(is.na(alpha))) stop("Value error: Missing values (NA/NaN) detected in matrix 'alpha'.")
  if (!is.null(m.temp) && any(is.na(m.temp))) stop("Value error: Missing values (NA/NaN) detected in 'm.temp'.")

  # Dispersal vector (d.spp) safeguards
  if (any(d.spp < 0)) stop("Value error: Dispersal coefficients in 'd.spp' cannot be negative.")
  if (sum(d.spp) == 0) stop("Mathematical error: 'd.spp' cannot sum to zero. Use a vector of 1s for a neutral effect or set 'd.spp = NULL'.")

  # Regional pool safeguards
  if (any(Meta.pool < 0)) stop("Value error: Relative abundances in 'Meta.pool' cannot be negative.")
  if (sum(Meta.pool) == 0) stop("Mathematical error: 'Meta.pool' cannot sum to zero.")

  # comm.fixed zero-sum crash protections
  if (!is.null(comm.fixed)) {
    if (any(comm.fixed < 0)) stop("Value error: Relative abundances in 'comm.fixed' cannot be negative.")
    if (is.matrix(comm.fixed)) {
      if (any(colSums(comm.fixed) == 0)) stop("Mathematical error: A column in 'comm.fixed' sums to zero.")
    } else {
      if (sum(comm.fixed) == 0) stop("Mathematical error: 'comm.fixed' cannot sum to zero.")
    }
  }

  # Ecological and Mathematical boundary safeguards
  if (m.pool < 0 || m.pool > 1)                     stop("'m.pool' regional immigration rate must be between 0 and 1.")
  if (prop.dead.by.it <= 0 || prop.dead.by.it >= 1) stop("'prop.dead.by.it' baseline mortality fraction must be between 0 and 1.")
  if (any(Js <= 0))                                 stop("Carrying capacities in 'Js' must be strictly positive integers.")
  if (nIterations <= 0)                             stop("Number of lottery iterations 'it' must be a positive integer.")

  # Enforce strict probability bounds on FF and catch all-zero crash conditions
  if (!is.null(FF)) {
    if (any(FF < 0) || any(FF > 1)) stop("Value error: All filtering coefficients in matrix 'FF' must scale between 0 and 1.")
    if (any(colSums(FF) == 0)) {
      stop("Mathematical error: 'FF' cannot contain a column/community of all zeros. This completely blocks recruitment and will crash the simulation. Use 1s for no effect.")
    }
  }

  # Competition matrix (alpha) safeguards
  if (!is.null(alpha)) {
    if (any(alpha < 0)) stop("Value error: Competition coefficients in 'alpha' cannot be negative.")
    if (all(alpha == 0)) {
      stop("Mathematical error: 'alpha' cannot be a matrix of all zeros (causes a 0/0 NaN loop crash). Set 'alpha = NULL' or use a matrix of 1s to remove the effect.")
    }
  }

  # Enforce Kelvin temperature protection
  if (any(Ts <= 0)) stop("Value error: Temperatures in 'Ts' must be strictly positive values expressed in Kelvin.")

}

