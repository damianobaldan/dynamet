#' Validate the inputs for masterEqMetacomm
#'
#' @inheritParams masterEqMetacomm
#' @keywords internal
#'
validateMetaInputs <- function(
    Meta.pool,
    d.spp,
    FF = NULL,
    m.pool,
    Js,
    M.migra = NULL,
    id.fixed = NULL,
    comm.fixed = NULL,
    init.comm,
    lottery,
    nIterations,
    prop.dead.by.it,
    id.obs = NULL,
    Ea,
    Ts,
    m.temp = NULL,
    alpha = NULL) {

  # ----------------------------------------------------------------------------
  # 1. CLASS AND TYPE VALIDATIONS
  # ----------------------------------------------------------------------------
  if (!is.numeric(m.pool) || length(m.pool) != 1) stop("'m.pool' must be a single numeric value.")
  if (!is.logical(lottery) || length(lottery) != 1) stop("'lottery' must be a single logical value.")
  if (!is.numeric(it) || length(it) != 1) stop("'it' must be a single numeric integer.")
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
  if (!is.null(id.obs) && !is.numeric(id.obs))     stop("'id.obs' contains invalid non-numeric values.")
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
  if (length(Ts) != C) stop(sprintf("Dimension mismatch: Temperature vector 'Ts' length (%d) must match 'Js' (%d).", length(Ts), C))

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
  if (!is.null(id.obs) && (any(id.obs < 1) || any(id.obs > C)))     stop("'id.obs' contains out-of-bounds community indices.")

  # ----------------------------------------------------------------------------
  # 3. VALUE AND BOUNDARY CONSTRAINTS (Includes Zero-Sum / NaN Protections)
  # ----------------------------------------------------------------------------
  # Missing data sweep
  all_inputs <- list(Meta.pool, d.spp, m.pool, Js, M.migra, it, prop.dead.by.it, Ea, Ts)
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
  if (it <= 0)                                      stop("Number of lottery iterations 'it' must be a positive integer.")

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
