#' Simulate Mutualistic Metacommunity Dynamics (Animals and Plants)
#'
#' Models bipartite metacommunity assembly and demographic turnover across a
#' network of local communities for two mutually interacting species groups
#' (animals and plants). Incorporates spatial migration, regional immigration,
#' species-specific dispersal, environmental filtering, intraspecific/interspecific
#' competition, temperature-dependent mortality, and bipartite facilitation.
#'
#' @param Meta.pool.animal Numeric vector of relative abundances for animal species (S_animal).
#' @param Meta.pool.plant Numeric vector of relative abundances for plant species (S_plant).
#' @param Js.animal Numeric vector of length C setting local carrying capacities for animals.
#' @param Js.plant Numeric vector of length C setting local carrying capacities for plants.
#' @param M.migra.animal Square matrix (C x C) of animal spatial migration probabilities.
#' @param M.migra.plant Square matrix (C x C) of plant spatial migration probabilities.
#' @param m.pool.animal Fraction (0 to 1) of animal recruits originating from regional pool.
#' @param m.pool.plant Fraction (0 to 1) of plant recruits originating from regional pool.
#' @param d.spp.animal Numeric vector (S_animal) of animal dispersal capabilities. Default is equal (1s).
#' @param d.spp.plant Numeric vector (S_plant) of plant dispersal capabilities. Default is equal (1s).
#' @param FF.animal Matrix (S_animal x C) of environmental filtering coefficients for animals.
#' @param FF.plant Matrix (S_plant x C) of environmental filtering coefficients for plants.
#' @param alpha.animal Square matrix (S_animal x S_animal) of animal competition coefficients.
#' @param alpha.plant Square matrix (S_plant x S_plant) of plant competition coefficients.
#' @param alpha.animal.plant Bipartite matrix (S_animal x S_plant) defining positive facilitation coefficients.
#' @param init.comm.animal Optional starting abundance matrix (S_animal x C) for animals.
#' @param init.comm.plant Optional starting abundance matrix (S_plant x C) for plants.
#' @param id.fixed Vector of community indices with static/locked compositions.
#' @param comm.fixed.animal Vector or matrix (S_animal x length(id.fixed)) of fixed animal profiles.
#' @param comm.fixed.plant Vector or matrix (S_plant x length(id.fixed)) of fixed plant profiles.
#' @param m.temp.animal Matrix or vector setting memory/seed bank weights for animals.
#' @param m.temp.plant Matrix or vector setting memory/seed bank weights for plants.
#' @param prop.dead.by.it Baseline turnover fraction for the coldest patch.
#' @param Ea Activation energy (eV) for Arrhenius temperature scaling.
#' @param Ts Vector of patch temperatures in Kelvin (> 0) or a single scalar temperature value.
#' @param lottery Logical. If \code{TRUE}, runs the iterative turnover phase.
#' @param nIterations Total timeline steps for the lottery loop.
#' @param verbose Logical. If \code{TRUE}, prints progress updates to console.
#'
#' @return A list with two matrices:
#' \item{Meta.animal}{Matrix (S_animal x C) of final animal abundances.}
#' \item{Meta.plant}{Matrix (S_plant x C) of final plant abundances.}
#'
#' @export
#'
#' @author Matias Arim & Ana Borthagaray
#'
masterEqMutualistic <- function(Meta.pool.animal, Meta.pool.plant,
                                Js.animal, Js.plant,
                                M.migra.animal, M.migra.plant,
                                m.pool.animal, m.pool.plant,
                                d.spp.animal = NULL, d.spp.plant = NULL,
                                FF.animal = NULL, FF.plant = NULL,
                                alpha.animal = NULL, alpha.plant = NULL, alpha.animal.plant = NULL,
                                init.comm.animal = NULL, init.comm.plant = NULL,
                                id.fixed = NULL, comm.fixed.animal = NULL, comm.fixed.plant = NULL,
                                m.temp.animal = 0, m.temp.plant = 0,
                                prop.dead.by.it = 0.05, Ea = 1e-05, Ts = 293.15,
                                lottery = TRUE, nIterations = 100, verbose = TRUE) {

  #----------- 1. VALIDATE INPUTS -----------
  validateMutualisticMetaInputs(
    Meta.pool.animal        = Meta.pool.animal,
    Meta.pool.plant         = Meta.pool.plant,
    Js.animal               = Js.animal,
    Js.plant                = Js.plant,
    M.migra.animal          = M.migra.animal,
    M.migra.plant           = M.migra.plant,
    m.pool.animal           = m.pool.animal,
    m.pool.plant            = m.pool.plant,
    d.spp.animal            = d.spp.animal,
    d.spp.plant             = d.spp.plant,
    FF.animal               = FF.animal,
    FF.plant                = FF.plant,
    alpha.animal            = alpha.animal,
    alpha.plant             = alpha.plant,
    alpha.animal.plant      = alpha.animal.plant,
    init.comm.animal        = init.comm.animal,
    init.comm.plant         = init.comm.plant,
    id.fixed                = id.fixed,
    comm.fixed.animal       = comm.fixed.animal,
    comm.fixed.plant        = comm.fixed.plant,
    m.temp.animal           = m.temp.animal,
    m.temp.plant            = m.temp.plant,
    prop.dead.by.it         = prop.dead.by.it,
    Ea                      = Ea,
    Ts                      = Ts,
    lottery                 = lottery,
    nIterations             = nIterations
  )

  #----------- 2. INITIALIZATION AND DATA NORMALIZATION -----------
  S.animal <- length(Meta.pool.animal)
  S.plant  <- length(Meta.pool.plant)
  C        <- length(Js.animal)

  # Dispersal vectors
  if (is.null(d.spp.animal)) d.spp.animal <- rep(1, S.animal)
  if (is.null(d.spp.plant))  d.spp.plant  <- rep(1, S.plant)
  d.spp.animal <- d.spp.animal / sum(d.spp.animal)
  d.spp.plant  <- d.spp.plant / sum(d.spp.plant)

  # Regional Pools
  Meta.pool.animal <- Meta.pool.animal / sum(Meta.pool.animal)
  Meta.pool.plant  <- Meta.pool.plant / sum(Meta.pool.plant)

  # Normalize Migration Matrices (Column sums normalized to 1)
  M.migra.animal <- sweep(M.migra.animal, 2, colSums(M.migra.animal), FUN = "/")
  M.migra.plant  <- sweep(M.migra.plant, 2, colSums(M.migra.plant), FUN = "/")

  # Fixed Communities Normalization
  if (!is.null(comm.fixed.animal)) {
    if (is.matrix(comm.fixed.animal)) {
      comm.fixed.animal <- sweep(comm.fixed.animal, 2, colSums(comm.fixed.animal), FUN = "/")
    } else {
      comm.fixed.animal <- comm.fixed.animal / sum(comm.fixed.animal)
    }
  }
  if (!is.null(comm.fixed.plant)) {
    if (is.matrix(comm.fixed.plant)) {
      comm.fixed.plant <- sweep(comm.fixed.plant, 2, colSums(comm.fixed.plant), FUN = "/")
    } else {
      comm.fixed.plant <- comm.fixed.plant / sum(comm.fixed.plant)
    }
  }

  # Environmental Filters
  if (is.null(FF.animal)) FF.animal <- matrix(1, nrow = S.animal, ncol = C)
  if (is.null(FF.plant))  FF.plant  <- matrix(1, nrow = S.plant, ncol = C)

  # Memory Matrices
  m.temp.animal <- formatMemoryMatrix(m.temp.animal, S.animal, C)
  m.temp.plant  <- formatMemoryMatrix(m.temp.plant, S.plant, C)

  # Expand Ts if input is scalar
  if (length(Ts) == 1) { Ts <- rep(Ts, C) }

  #----------- 3. THERMAL DEPENDENCE CALCULATIONS -----------
  min.dead.Tmin     <- prop.dead.by.it / exp(-Ea / (min(Ts) * 8.62e-05))
  prop.dead.by.comm <- min.dead.Tmin * exp(-Ea / (Ts * 8.62e-05))

  if (max(prop.dead.by.comm) > 0.95) {
    stop("Error: Temperature-driven mortality replaces >95% of individuals. Decrease prop.dead.by.it.")
  }

  dead.by.it.animal <- round(prop.dead.by.comm * Js.animal, 0)
  dead.by.it.animal <- ifelse(dead.by.it.animal < 2, 2, dead.by.it.animal)

  dead.by.it.plant  <- round(prop.dead.by.comm * Js.plant, 0)
  dead.by.it.plant  <- ifelse(dead.by.it.plant < 2, 2, dead.by.it.plant)

  if (!is.null(id.fixed)) {
    dead.by.it.animal[id.fixed] <- 0
    dead.by.it.plant[id.fixed]  <- 0
  }

  max.dead.by.it.animal <- max(dead.by.it.animal)
  max.dead.by.it.plant  <- max(dead.by.it.plant)

  #----------- 4. COALESCENT ASSEMBLY PHASE -----------
  if (is.null(init.comm.animal)) {
    Meta.animal <- matrix(0, nrow = S.animal, ncol = C)
    for (i in 1:C) {
      Meta.animal[, i] <- stats::rmultinom(1, 1, Meta.pool.animal * d.spp.animal * FF.animal[, i])
    }
  } else {
    Meta.animal <- init.comm.animal
  }

  if (is.null(init.comm.plant)) {
    Meta.plant <- matrix(0, nrow = S.plant, ncol = C)
    for (i in 1:C) {
      Meta.plant[, i] <- stats::rmultinom(1, 1, Meta.pool.plant * d.spp.plant * FF.plant[, i])
    }
  } else {
    Meta.plant <- init.comm.plant
  }

  max_capacity <- max(c(Js.animal, Js.plant))

  for (ii in 2:max_capacity) {

    id.j.animal <- which(colSums(Meta.animal) < Js.animal & colSums(Meta.animal) < ii)
    id.j.plant  <- which(colSums(Meta.plant) < Js.plant & colSums(Meta.plant) < ii)

    if (!is.null(id.fixed)) {
      id.j.animal <- setdiff(id.j.animal, id.fixed)
      id.j.plant  <- setdiff(id.j.plant, id.fixed)
    }

    # Lock fixed communities at proportional capacity
    if (!is.null(id.fixed)) {
      for (idx in seq_along(id.fixed)) {
        f_id <- id.fixed[idx]
        prof.a <- if (is.matrix(comm.fixed.animal)) comm.fixed.animal[, idx] else comm.fixed.animal
        prof.p <- if (is.matrix(comm.fixed.plant))  comm.fixed.plant[, idx]  else comm.fixed.plant
        Meta.animal[, f_id] <- prof.a * min(ii - 1, Js.animal[f_id])
        Meta.plant[, f_id]  <- prof.p * min(ii - 1, Js.plant[f_id])
      }
    }

    # --- ANIMAL RECRUITMENT ---
    if (length(id.j.animal) > 0) {
      Pool.neighbor.animal <- (Meta.animal %*% M.migra.animal) * d.spp.animal * FF.animal

      if (!is.null(alpha.animal)) {
        overlap.animal <- alpha.animal %*% Meta.animal
        col_sums_overlap <- colSums(overlap.animal)
        col_sums_overlap[col_sums_overlap == 0] <- 1
        overlap.animal <- sweep(overlap.animal, 2, col_sums_overlap, FUN = "/")
        Pool.neighbor.animal <- Pool.neighbor.animal * (1 - overlap.animal)
      }

      if (!is.null(alpha.animal.plant)) {
        facil.animal <- (alpha.animal.plant %*% Meta.plant) + 1
        col_sums_facil <- colSums(facil.animal)
        col_sums_facil[col_sums_facil == 0] <- 1
        facil.animal <- sweep(facil.animal, 2, col_sums_facil, FUN = "/")
        Pool.neighbor.animal <- Pool.neighbor.animal * facil.animal
      }

      zero_cols <- which(colSums(Pool.neighbor.animal) == 0)
      for (zc in zero_cols) Pool.neighbor.animal[, zc] <- Meta.pool.animal

      if (length(id.j.animal) > 1) {
        new.a <- apply(Pool.neighbor.animal[, id.j.animal, drop = FALSE], 2, born, dead.by.it = 1, M.pool = Meta.pool.animal, m.pool = m.pool.animal)
        Meta.animal[, id.j.animal] <- Meta.animal[, id.j.animal] + new.a
      } else {
        Meta.animal[, id.j.animal] <- Meta.animal[, id.j.animal] + born(probs = Pool.neighbor.animal[, id.j.animal], dead.by.it = 1, M.pool = Meta.pool.animal, m.pool = m.pool.animal)
      }
    }

    # --- PLANT RECRUITMENT ---
    if (length(id.j.plant) > 0) {
      Pool.neighbor.plant <- (Meta.plant %*% M.migra.plant) * d.spp.plant * FF.plant

      if (!is.null(alpha.plant)) {
        overlap.plant <- alpha.plant %*% Meta.plant
        col_sums_overlap <- colSums(overlap.plant)
        col_sums_overlap[col_sums_overlap == 0] <- 1
        overlap.plant <- sweep(overlap.plant, 2, col_sums_overlap, FUN = "/")
        Pool.neighbor.plant <- Pool.neighbor.plant * (1 - overlap.plant)
      }

      if (!is.null(alpha.animal.plant)) {
        facil.plant <- (t(alpha.animal.plant) %*% Meta.animal) + 1
        col_sums_facil <- colSums(facil.plant)
        col_sums_facil[col_sums_facil == 0] <- 1
        facil.plant <- sweep(facil.plant, 2, col_sums_facil, FUN = "/")
        Pool.neighbor.plant <- Pool.neighbor.plant * facil.plant
      }

      zero_cols <- which(colSums(Pool.neighbor.plant) == 0)
      for (zc in zero_cols) Pool.neighbor.plant[, zc] <- Meta.pool.plant

      if (length(id.j.plant) > 1) {
        new.p <- apply(Pool.neighbor.plant[, id.j.plant, drop = FALSE], 2, born, dead.by.it = 1, M.pool = Meta.pool.plant, m.pool = m.pool.plant)
        Meta.plant[, id.j.plant] <- Meta.plant[, id.j.plant] + new.p
      } else {
        Meta.plant[, id.j.plant] <- Meta.plant[, id.j.plant] + born(probs = Pool.neighbor.plant[, id.j.plant], dead.by.it = 1, M.pool = Meta.pool.plant, m.pool = m.pool.plant)
      }
    }

    if (verbose) cat("Coalescent construction step J:", ii, "of", max_capacity, "\n")
  }

  #----------- 5. LOTTERY DYNAMICS PHASE -----------
  if (lottery) {

    Meta.lag.animal <- Meta.animal
    Meta.lag.plant  <- Meta.plant

    if (!is.null(id.fixed)) {
      for (idx in seq_along(id.fixed)) {
        f_id <- id.fixed[idx]
        prof.a <- if (is.matrix(comm.fixed.animal)) comm.fixed.animal[, idx] else comm.fixed.animal
        prof.p <- if (is.matrix(comm.fixed.plant))  comm.fixed.plant[, idx]  else comm.fixed.plant
        Meta.animal[, f_id] <- round(prof.a * Js.animal[f_id], 0)
        Meta.plant[, f_id]  <- round(prof.p * Js.plant[f_id], 0)
      }
    }

    generations <- seq(1, nIterations, 1 / prop.dead.by.it)

    for (iteration in 1:nIterations) {

      if (iteration %in% generations) {
        Meta.lag.animal <- Meta.animal
        Meta.lag.plant  <- Meta.plant
      }

      # --- Death Sub-phase: Animals ---
      for (dead in 1:max.dead.by.it.animal) {
        id.dead <- which(dead.by.it.animal >= dead)
        if (length(id.dead) > 1) {
          Meta.animal[, id.dead] <- Meta.animal[, id.dead] - apply(Meta.animal[, id.dead, drop = FALSE] * (1.001 - FF.animal[, id.dead, drop = FALSE]), 2, FUN = change, change = 1)
        } else if (length(id.dead) == 1) {
          Meta.animal[, id.dead] <- Meta.animal[, id.dead] - change(probs = Meta.animal[, id.dead] * (1.001 - FF.animal[, id.dead]), change = 1)
        }
      }

      # --- Death Sub-phase: Plants ---
      for (dead in 1:max.dead.by.it.plant) {
        id.dead <- which(dead.by.it.plant >= dead)
        if (length(id.dead) > 1) {
          Meta.plant[, id.dead] <- Meta.plant[, id.dead] - apply(Meta.plant[, id.dead, drop = FALSE] * (1.001 - FF.plant[, id.dead, drop = FALSE]), 2, FUN = change, change = 1)
        } else if (length(id.dead) == 1) {
          Meta.plant[, id.dead] <- Meta.plant[, id.dead] - change(probs = Meta.plant[, id.dead] * (1.001 - FF.plant[, id.dead]), change = 1)
        }
      }

      # --- Recruitment Sub-phase: Animals ---
      Pool.neighbor.animal <- (Meta.animal %*% M.migra.animal) * d.spp.animal * FF.animal
      Pool.neighbor.animal <- Pool.neighbor.animal * (1 - m.temp.animal) + Meta.lag.animal * m.temp.animal

      if (!is.null(alpha.animal)) {
        overlap.animal <- alpha.animal %*% Meta.animal
        col_sums_overlap <- colSums(overlap.animal)
        col_sums_overlap[col_sums_overlap == 0] <- 1
        overlap.animal <- sweep(overlap.animal, 2, col_sums_overlap, FUN = "/")
        Pool.neighbor.animal <- Pool.neighbor.animal * (1 - overlap.animal)
      }

      if (!is.null(alpha.animal.plant)) {
        facil.animal <- (alpha.animal.plant %*% Meta.plant) + 1
        col_sums_facil <- colSums(facil.animal)
        col_sums_facil[col_sums_facil == 0] <- 1
        facil.animal <- sweep(facil.animal, 2, col_sums_facil, FUN = "/")
        Pool.neighbor.animal <- Pool.neighbor.animal * facil.animal
      }

      col_sums.a <- colSums(Pool.neighbor.animal)
      col_sums.a[col_sums.a == 0] <- 1
      Pool.norm.animal <- sweep(Pool.neighbor.animal, 2, col_sums.a, FUN = "/")
      Prob.mat.animal  <- (1 - m.pool.animal) * Pool.norm.animal + m.pool.animal * Meta.pool.animal

      for (i in seq_along(dead.by.it.animal)) {
        if (dead.by.it.animal[i] > 0) {
          Meta.animal[, i] <- Meta.animal[, i] + stats::rmultinom(1, size = dead.by.it.animal[i], prob = Prob.mat.animal[, i])
        }
      }

      # --- Recruitment Sub-phase: Plants ---
      Pool.neighbor.plant <- (Meta.plant %*% M.migra.plant) * d.spp.plant * FF.plant
      Pool.neighbor.plant <- Pool.neighbor.plant * (1 - m.temp.plant) + Meta.lag.plant * m.temp.plant

      if (!is.null(alpha.plant)) {
        overlap.plant <- alpha.plant %*% Meta.plant
        col_sums_overlap <- colSums(overlap.plant)
        col_sums_overlap[col_sums_overlap == 0] <- 1
        overlap.plant <- sweep(overlap.plant, 2, col_sums_overlap, FUN = "/")
        Pool.neighbor.plant <- Pool.neighbor.plant * (1 - overlap.plant)
      }

      if (!is.null(alpha.animal.plant)) {
        facil.plant <- (t(alpha.animal.plant) %*% Meta.animal) + 1
        col_sums_facil <- colSums(facil.plant)
        col_sums_facil[col_sums_facil == 0] <- 1
        facil.plant <- sweep(facil.plant, 2, col_sums_facil, FUN = "/")
        Pool.neighbor.plant <- Pool.neighbor.plant * facil.plant
      }

      col_sums.p <- colSums(Pool.neighbor.plant)
      col_sums.p[col_sums.p == 0] <- 1
      Pool.norm.plant <- sweep(Pool.neighbor.plant, 2, col_sums.p, FUN = "/")
      Prob.mat.plant  <- (1 - m.pool.plant) * Pool.norm.plant + m.pool.plant * Meta.pool.plant

      for (i in seq_along(dead.by.it.plant)) {
        if (dead.by.it.plant[i] > 0) {
          Meta.plant[, i] <- Meta.plant[, i] + stats::rmultinom(1, size = dead.by.it.plant[i], prob = Prob.mat.plant[, i])
        }
      }

      if (verbose) cat("Lottery iteration", iteration, "of", nIterations, "\n")
    }
  }

  return(list(Meta.animal = Meta.animal, Meta.plant = Meta.plant))
}




#' Validate the inputs for masterEqMetacommMutualistic
#'
#' @inheritParams masterEqMutualistic
#' @keywords internal
#'
validateMutualisticMetaInputs <- function(Meta.pool.animal,
                                          Meta.pool.plant,
                                          Js.animal,
                                          Js.plant,
                                          M.migra.animal,
                                          M.migra.plant,
                                          m.pool.animal,
                                          m.pool.plant,
                                          d.spp.animal = NULL,
                                          d.spp.plant = NULL,
                                          FF.animal = NULL,
                                          FF.plant = NULL,
                                          alpha.animal = NULL,
                                          alpha.plant = NULL,
                                          alpha.animal.plant = NULL,
                                          init.comm.animal = NULL,
                                          init.comm.plant = NULL,
                                          id.fixed = NULL,
                                          comm.fixed.animal = NULL,
                                          comm.fixed.plant = NULL,
                                          m.temp.animal = NULL,
                                          m.temp.plant = NULL,
                                          prop.dead.by.it,
                                          Ea,
                                          Ts,
                                          lottery,
                                          nIterations) {

  # ----------------------------------------------------------------------------
  # 1. CLASS AND TYPE VALIDATIONS
  # ----------------------------------------------------------------------------
  if (!is.numeric(m.pool.animal) || length(m.pool.animal) != 1) stop("'m.pool.animal' must be a single numeric value.")
  if (!is.numeric(m.pool.plant) || length(m.pool.plant) != 1)   stop("'m.pool.plant' must be a single numeric value.")
  if (!is.logical(lottery) || length(lottery) != 1)            stop("'lottery' must be a single logical value.")
  if (!is.numeric(nIterations) || length(nIterations) != 1)    stop("'nIterations' must be a single numeric integer.")
  if (!is.numeric(prop.dead.by.it) || length(prop.dead.by.it) != 1) stop("'prop.dead.by.it' must be a single numeric fraction.")
  if (!is.numeric(Ea) || length(Ea) != 1)                      stop("'Ea' must be a single numeric value.")

  if (!is.numeric(Meta.pool.animal)) stop("'Meta.pool.animal' must be a numeric vector.")
  if (!is.numeric(Meta.pool.plant))  stop("'Meta.pool.plant' must be a numeric vector.")
  if (!is.numeric(Js.animal))        stop("'Js.animal' must be a numeric vector.")
  if (!is.numeric(Js.plant))         stop("'Js.plant' must be a numeric vector.")
  if (!is.numeric(Ts))               stop("'Ts' must be a numeric vector.")

  if (!is.null(d.spp.animal) && !is.numeric(d.spp.animal)) stop("'d.spp.animal' must be a numeric vector.")
  if (!is.null(d.spp.plant) && !is.numeric(d.spp.plant))   stop("'d.spp.plant' must be a numeric vector.")

  if (!is.null(FF.animal) && (!is.matrix(FF.animal) || !is.numeric(FF.animal))) stop("'FF.animal' must be a numeric matrix.")
  if (!is.null(FF.plant) && (!is.matrix(FF.plant) || !is.numeric(FF.plant)))   stop("'FF.plant' must be a numeric matrix.")

  if (!is.null(alpha.animal) && (!is.matrix(alpha.animal) || !is.numeric(alpha.animal)))             stop("'alpha.animal' must be a numeric matrix.")
  if (!is.null(alpha.plant) && (!is.matrix(alpha.plant) || !is.numeric(alpha.plant)))               stop("'alpha.plant' must be a numeric matrix.")
  if (!is.null(alpha.animal.plant) && (!is.matrix(alpha.animal.plant) || !is.numeric(alpha.animal.plant))) stop("'alpha.animal.plant' must be a numeric matrix.")

  if (!is.null(init.comm.animal) && (!is.matrix(init.comm.animal) || !is.numeric(init.comm.animal))) stop("'init.comm.animal' must be a numeric matrix.")
  if (!is.null(init.comm.plant) && (!is.matrix(init.comm.plant) || !is.numeric(init.comm.plant)))   stop("'init.comm.plant' must be a numeric matrix.")

  if (!is.null(M.migra.animal) && (!is.matrix(M.migra.animal) || !is.numeric(M.migra.animal))) stop("'M.migra.animal' must be a numeric matrix.")
  if (!is.null(M.migra.plant) && (!is.matrix(M.migra.plant) || !is.numeric(M.migra.plant)))   stop("'M.migra.plant' must be a numeric matrix.")

  # Configuration cross-dependence safety
  if (!is.null(id.fixed) && (is.null(comm.fixed.animal) || is.null(comm.fixed.plant))) {
    stop("Configuration error: 'id.fixed' was provided, but 'comm.fixed.animal' or 'comm.fixed.plant' is NULL.")
  }
  if (is.null(id.fixed) && (!is.null(comm.fixed.animal) || !is.null(comm.fixed.plant))) {
    stop("Configuration error: Fixed profiles were provided, but 'id.fixed' is NULL.")
  }

  if (!is.null(comm.fixed.animal) && !is.numeric(comm.fixed.animal)) stop("'comm.fixed.animal' must be a numeric vector or matrix.")
  if (!is.null(comm.fixed.plant) && !is.numeric(comm.fixed.plant))   stop("'comm.fixed.plant' must be a numeric vector or matrix.")

  if (!is.null(id.fixed) && !is.numeric(id.fixed)) stop("'id.fixed' must be a numeric vector.")
  if (!is.null(m.temp.animal) && !is.numeric(m.temp.animal)) stop("'m.temp.animal' must be a numeric vector or matrix.")
  if (!is.null(m.temp.plant) && !is.numeric(m.temp.plant))   stop("'m.temp.plant' must be a numeric vector or matrix.")

  # ----------------------------------------------------------------------------
  # 2. DIMENSION CONSISTENCY VALIDATIONS
  # ----------------------------------------------------------------------------
  S.a <- length(Meta.pool.animal)
  S.p <- length(Meta.pool.plant)
  C.a <- length(Js.animal)
  C.p <- length(Js.plant)

  if (C.a != C.p) stop(sprintf("Dimension mismatch: 'Js.animal' length (%d) must match 'Js.plant' (%d).", C.a, C.p))
  C <- C.a

  if (!(length(Ts) %in% c(1, C))) {
    stop(sprintf("Dimension mismatch: 'Ts' length (%d) must match total communities C (%d).", length(Ts), C))
  }

  # Species-dimension alignments (Animals)
  if (!is.null(d.spp.animal) && length(d.spp.animal) != S.a) {
    stop(sprintf("Dimension mismatch: 'd.spp.animal' length (%d) must match 'Meta.pool.animal' (%d).", length(d.spp.animal), S.a))
  }
  if (!is.null(FF.animal) && nrow(FF.animal) != S.a) {
    stop(sprintf("Dimension mismatch: Environmental filter matrix 'FF.animal' must have %d rows (Animal Species).", S.a))
  }
  if (!is.null(alpha.animal) && (nrow(alpha.animal) != S.a || ncol(alpha.animal) != S.a)) {
    stop(sprintf("Dimension mismatch: Interspecific competition matrix 'alpha.animal' must be a square matrix of %d x %d.", S.a, S.a))
  }
  if (!is.null(init.comm.animal) && (nrow(init.comm.animal) != S.a || ncol(init.comm.animal) != C)) {
    stop(sprintf("Dimension mismatch: Matrix 'init.comm.animal' must be of dimension %d x %d.", S.a, C))
  }

  # Species-dimension alignments (Plants)
  if (!is.null(d.spp.plant) && length(d.spp.plant) != S.p) {
    stop(sprintf("Dimension mismatch: 'd.spp.plant' length (%d) must match 'Meta.pool.plant' (%d).", length(d.spp.plant), S.p))
  }
  if (!is.null(FF.plant) && nrow(FF.plant) != S.p) {
    stop(sprintf("Dimension mismatch: Environmental filter matrix 'FF.plant' must have %d rows (Plant Species).", S.p))
  }
  if (!is.null(alpha.plant) && (nrow(alpha.plant) != S.p || ncol(alpha.plant) != S.p)) {
    stop(sprintf("Dimension mismatch: Interspecific competition matrix 'alpha.plant' must be a square matrix of %d x %d.", S.p, S.p))
  }
  if (!is.null(init.comm.plant) && (nrow(init.comm.plant) != S.p || ncol(init.comm.plant) != C)) {
    stop(sprintf("Dimension mismatch: Matrix 'init.comm.plant' must be of dimension %d x %d.", S.p, C))
  }

  # Bipartite Mutualism Interaction Matrix Alignment
  if (!is.null(alpha.animal.plant) && (nrow(alpha.animal.plant) != S.a || ncol(alpha.animal.plant) != S.p)) {
    stop(sprintf("Dimension mismatch: Mutualism matrix 'alpha.animal.plant' must have dimensions S_animal (%d) x S_plant (%d).", S.a, S.p))
  }

  # Community-dimension alignments (Columns)
  if (!is.null(FF.animal) && ncol(FF.animal) != C) stop(sprintf("Dimension mismatch: Environmental filter matrix 'FF.animal' must have %d columns (Communities).", C))
  if (!is.null(FF.plant) && ncol(FF.plant) != C)   stop(sprintf("Dimension mismatch: Environmental filter matrix 'FF.plant' must have %d columns (Communities).", C))

  if (!is.null(M.migra.animal)) {
    if (nrow(M.migra.animal) != C || ncol(M.migra.animal) != C) stop(sprintf("Dimension mismatch: Migration matrix 'M.migra.animal' must be a square matrix of %d x %d.", C, C))
  } else {
    stop("Missing required input: 'M.migra.animal' matrix cannot be NULL.")
  }

  if (!is.null(M.migra.plant)) {
    if (nrow(M.migra.plant) != C || ncol(M.migra.plant) != C) stop(sprintf("Dimension mismatch: Migration matrix 'M.migra.plant' must be a square matrix of %d x %d.", C, C))
  } else {
    stop("Missing required input: 'M.migra.plant' matrix cannot be NULL.")
  }

  # Fixed profiles dimensions
  if (!is.null(comm.fixed.animal)) {
    if (is.matrix(comm.fixed.animal)) {
      if (nrow(comm.fixed.animal) != S.a || ncol(comm.fixed.animal) != length(id.fixed)) {
        stop(sprintf("Dimension mismatch: 'comm.fixed.animal' matrix must have dimensions S_animal (%d) x length(id.fixed) (%d).", S.a, length(id.fixed)))
      }
    } else {
      if (length(comm.fixed.animal) != S.a) {
        stop(sprintf("Dimension mismatch: 'comm.fixed.animal' vector must be of length S_animal (%d).", S.a))
      }
    }
  }

  if (!is.null(comm.fixed.plant)) {
    if (is.matrix(comm.fixed.plant)) {
      if (nrow(comm.fixed.plant) != S.p || ncol(comm.fixed.plant) != length(id.fixed)) {
        stop(sprintf("Dimension mismatch: 'comm.fixed.plant' matrix must have dimensions S_plant (%d) x length(id.fixed) (%d).", S.p, length(id.fixed)))
      }
    } else {
      if (length(comm.fixed.plant) != S.p) {
        stop(sprintf("Dimension mismatch: 'comm.fixed.plant' vector must be of length S_plant (%d).", S.p))
      }
    }
  }

  # Memory Matrix/Vector dimensions
  if (!is.null(m.temp.animal)) {
    if (is.matrix(m.temp.animal)) {
      if (nrow(m.temp.animal) != S.a || ncol(m.temp.animal) != C) {
        stop(sprintf("Dimension mismatch: Matrix 'm.temp.animal' must match landscape dimensions (%d rows x %d columns).", S.a, C))
      }
    } else if (length(m.temp.animal) > 1 && length(m.temp.animal) != C) {
      stop(sprintf("Dimension mismatch: Vector 'm.temp.animal' length (%d) must match the number of communities (%d).", length(m.temp.animal), C))
    }
  }

  if (!is.null(m.temp.plant)) {
    if (is.matrix(m.temp.plant)) {
      if (nrow(m.temp.plant) != S.p || ncol(m.temp.plant) != C) {
        stop(sprintf("Dimension mismatch: Matrix 'm.temp.plant' must match landscape dimensions (%d rows x %d columns).", S.p, C))
      }
    } else if (length(m.temp.plant) > 1 && length(m.temp.plant) != C) {
      stop(sprintf("Dimension mismatch: Vector 'm.temp.plant' length (%d) must match the number of communities (%d).", length(m.temp.plant), C))
    }
  }

  if (!is.null(id.fixed) && (any(id.fixed < 1) || any(id.fixed > C))) {
    stop("'id.fixed' contains out-of-bounds community indices.")
  }

  # ----------------------------------------------------------------------------
  # 3. VALUE AND BOUNDARY CONSTRAINTS (Includes Zero-Sum / NaN Protections)
  # ----------------------------------------------------------------------------
  # Missing data sweep across primary parameters
  all_inputs <- list(Meta.pool.animal, Meta.pool.plant, m.pool.animal, m.pool.plant,
                     Js.animal, Js.plant, M.migra.animal, M.migra.plant,
                     nIterations, prop.dead.by.it, Ea, Ts)
  if (any(sapply(all_inputs, function(x) any(is.na(x))))) stop("Value error: Missing values (NA/NaN) detected in primary inputs.")

  if (!is.null(d.spp.animal) && any(is.na(d.spp.animal)))           stop("Value error: Missing values (NA/NaN) detected in 'd.spp.animal'.")
  if (!is.null(d.spp.plant) && any(is.na(d.spp.plant)))             stop("Value error: Missing values (NA/NaN) detected in 'd.spp.plant'.")
  if (!is.null(FF.animal) && any(is.na(FF.animal)))                 stop("Value error: Missing values (NA/NaN) detected in matrix 'FF.animal'.")
  if (!is.null(FF.plant) && any(is.na(FF.plant)))                   stop("Value error: Missing values (NA/NaN) detected in matrix 'FF.plant'.")
  if (!is.null(alpha.animal) && any(is.na(alpha.animal)))           stop("Value error: Missing values (NA/NaN) detected in matrix 'alpha.animal'.")
  if (!is.null(alpha.plant) && any(is.na(alpha.plant)))             stop("Value error: Missing values (NA/NaN) detected in matrix 'alpha.plant'.")
  if (!is.null(alpha.animal.plant) && any(is.na(alpha.animal.plant))) stop("Value error: Missing values (NA/NaN) detected in matrix 'alpha.animal.plant'.")
  if (!is.null(m.temp.animal) && any(is.na(m.temp.animal)))         stop("Value error: Missing values (NA/NaN) detected in 'm.temp.animal'.")
  if (!is.null(m.temp.plant) && any(is.na(m.temp.plant)))           stop("Value error: Missing values (NA/NaN) detected in 'm.temp.plant'.")

  # Dispersal vector safeguards
  if (!is.null(d.spp.animal)) {
    if (any(d.spp.animal < 0)) stop("Value error: Dispersal coefficients in 'd.spp.animal' cannot be negative.")
    if (sum(d.spp.animal) == 0) stop("Mathematical error: 'd.spp.animal' cannot sum to zero. Set 'd.spp.animal = NULL' for equal dispersal.")
  }
  if (!is.null(d.spp.plant)) {
    if (any(d.spp.plant < 0)) stop("Value error: Dispersal coefficients in 'd.spp.plant' cannot be negative.")
    if (sum(d.spp.plant) == 0) stop("Mathematical error: 'd.spp.plant' cannot sum to zero. Set 'd.spp.plant = NULL' for equal dispersal.")
  }

  # Regional pool safeguards
  if (any(Meta.pool.animal < 0)) stop("Value error: Relative abundances in 'Meta.pool.animal' cannot be negative.")
  if (sum(Meta.pool.animal) == 0) stop("Mathematical error: 'Meta.pool.animal' cannot sum to zero.")
  if (any(Meta.pool.plant < 0))  stop("Value error: Relative abundances in 'Meta.pool.plant' cannot be negative.")
  if (sum(Meta.pool.plant) == 0)  stop("Mathematical error: 'Meta.pool.plant' cannot sum to zero.")

  # comm.fixed safeguards
  if (!is.null(comm.fixed.animal)) {
    if (any(comm.fixed.animal < 0)) stop("Value error: Relative abundances in 'comm.fixed.animal' cannot be negative.")
    if (is.matrix(comm.fixed.animal)) {
      if (any(colSums(comm.fixed.animal) == 0)) stop("Mathematical error: A column in 'comm.fixed.animal' sums to zero.")
    } else {
      if (sum(comm.fixed.animal) == 0) stop("Mathematical error: 'comm.fixed.animal' cannot sum to zero.")
    }
  }
  if (!is.null(comm.fixed.plant)) {
    if (any(comm.fixed.plant < 0)) stop("Value error: Relative abundances in 'comm.fixed.plant' cannot be negative.")
    if (is.matrix(comm.fixed.plant)) {
      if (any(colSums(comm.fixed.plant) == 0)) stop("Mathematical error: A column in 'comm.fixed.plant' sums to zero.")
    } else {
      if (sum(comm.fixed.plant) == 0) stop("Mathematical error: 'comm.fixed.plant' cannot sum to zero.")
    }
  }

  # Migration matrix safeguards
  if (any(M.migra.animal < 0)) stop("Value error: Migration probabilities in 'M.migra.animal' cannot be negative.")
  if (any(colSums(M.migra.animal) == 0)) stop("Mathematical error: 'M.migra.animal' contains a column with a sum of zero.")
  if (any(M.migra.plant < 0))  stop("Value error: Migration probabilities in 'M.migra.plant' cannot be negative.")
  if (any(colSums(M.migra.plant) == 0))  stop("Mathematical error: 'M.migra.plant' contains a column with a sum of zero.")

  # Ecological and Mathematical boundary safeguards
  if (m.pool.animal < 0 || m.pool.animal > 1)       stop("'m.pool.animal' regional immigration rate must be between 0 and 1.")
  if (m.pool.plant < 0 || m.pool.plant > 1)         stop("'m.pool.plant' regional immigration rate must be between 0 and 1.")
  if (prop.dead.by.it <= 0 || prop.dead.by.it >= 1) stop("'prop.dead.by.it' baseline mortality fraction must be strictly between 0 and 1.")
  if (any(Js.animal <= 0))                          stop("Carrying capacities in 'Js.animal' must be strictly positive integers.")
  if (any(Js.plant <= 0))                           stop("Carrying capacities in 'Js.plant' must be strictly positive integers.")
  if (nIterations <= 0)                             stop("Number of lottery iterations 'nIterations' must be a positive integer.")

  # Environmental filter bounds and zero-column protection
  if (!is.null(FF.animal)) {
    if (any(FF.animal < 0) || any(FF.animal > 1)) stop("Value error: Filtering coefficients in 'FF.animal' must scale between 0 and 1.")
    if (any(colSums(FF.animal) == 0)) {
      stop("Mathematical error: 'FF.animal' cannot contain a column of all zeros (completely blocks recruitment).")
    }
  }
  if (!is.null(FF.plant)) {
    if (any(FF.plant < 0) || any(FF.plant > 1)) stop("Value error: Filtering coefficients in 'FF.plant' must scale between 0 and 1.")
    if (any(colSums(FF.plant) == 0)) {
      stop("Mathematical error: 'FF.plant' cannot contain a column of all zeros (completely blocks recruitment).")
    }
  }

  # Competition matrix safeguards
  if (!is.null(alpha.animal)) {
    if (any(alpha.animal < 0)) stop("Value error: Competition coefficients in 'alpha.animal' cannot be negative.")
    if (all(alpha.animal == 0)) {
      stop("Mathematical error: 'alpha.animal' cannot be all zeros (causes a NaN loop crash). Set 'alpha.animal = NULL' to omit.")
    }
  }
  if (!is.null(alpha.plant)) {
    if (any(alpha.plant < 0)) stop("Value error: Competition coefficients in 'alpha.plant' cannot be negative.")
    if (all(alpha.plant == 0)) {
      stop("Mathematical error: 'alpha.plant' cannot be all zeros (causes a NaN loop crash). Set 'alpha.plant = NULL' to omit.")
    }
  }

  # Mutualism interaction safeguards
  if (!is.null(alpha.animal.plant)) {
    if (any(alpha.animal.plant < 0)) stop("Value error: Interaction coefficients in 'alpha.animal.plant' cannot be negative.")
  }

  # Temperature protection
  if (any(Ts <= 0)) stop("Value error: Temperatures in 'Ts' must be strictly positive values expressed in Kelvin.")

}




#' Format Memory Matrices
#' @keywords internal
formatMemoryMatrix <- function(m.temp, S, C) {
  if (is.null(m.temp)) {
    return(matrix(0, nrow = S, ncol = C))
  }
  if (is.matrix(m.temp)) {
    return(m.temp)
  }
  if (length(m.temp) == 1) {
    return(matrix(m.temp, nrow = S, ncol = C))
  }
  if (length(m.temp) == C) {
    return(matrix(rep(m.temp, each = S), nrow = S, ncol = C))
  }
  stop("Invalid dimensions for memory parameter 'm.temp'.")
}
