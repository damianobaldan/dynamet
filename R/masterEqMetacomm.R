#' Simulate Metacommunity Dynamics via Coalescent Assembly and lottery Phases
#'
#' @description
#' This function models metacommunity assembly and demographic turnover across a
#' network of local communities. It initiates patches using a coalescent
#' process up to local carrying capacities (\code{Js}). Then, an optional neutral/niche
#' lottery dynamic is used to simulate metacommunities. The simulation incorporates spatial migration, regional pool
#' immigration, self-recruitment, species-specific dispersal constraints, environmental filtering,
#' interspecific competition (\code{Alfa}), and temperature-dependent mortality
#' scaled via the Arrhenius equation.
#'
#' @param Meta.pool A numeric vector of length S representing the relative abundances
#'   or probabilities of species within the regional species pool.
#' @param d.spp A numeric vector of length S dictating species-specific dispersal traits
#'   or regional recruitment constraints.
#' @param FF An optional numeric matrix of dimensions S x C representing local environmental
#'   filtering filters. Individual coefficients must range between 0 and 1, where 1 indicates
#'   perfect environmental match (no filtering penalty). If \code{NULL} (default), no environmental filtering is applied.
#' @param m.pool A single numeric value between 0 and 1 defining the probability of
#'   recruitment originating from the global regional pool rather than local/neighboring sources.
#' @param Js A numeric vector of length C setting the local carrying capacity (total
#'   individual slots) for each community patch.
#' @param M.migra A square numeric matrix of dimensions C x C establishing spatial
#'   migration connectivity and dispersal probabilities between patches. Cannot be NULL.
#' @param id.fixed An optional numeric vector containing indices of communities whose
#'   compositions are static and locked during simulation steps.
#' @param comm.fixed An optional numeric vector of length S outlining the fixed relative
#'   abundance profile assigned to the patches listed in \code{id.fixed}. These relative
#'   abundances are automatically scaled to the specific local carrying capacity (\code{Js})
#'   of each fixed patch.
#' @param init.comm An optional numeric matrix of dimensions S x C representing the custom
#'   starting abundance counts for all species across patches. Required if \code{coalescence = FALSE}.
#' @param lottery A single logical value. If \code{TRUE}, triggers the demographic turnover
#'   iterative lottery phase after initial community assembly.
#' @param it A single numeric integer specifying the total timeline steps/iterations to
#'   run inside the lottery loop.
#' @param prop.dead.by.it A single numeric fraction (0, 1) defining the baseline mortality
#'   turnover rate, anchored at the coldest patch.
#' @param Ea A single numeric value representing the activation energy (in eV) utilized to
#'   scale temperature-driven mortality.
#' @param Ts A numeric vector of length C specifying the local patch ambient temperatures
#'   expressed strictly in Kelvin units (> 0).
#' @param m.temp An optional numeric vector of length C, or a matrix of dimensions S x C,
#'   dictating community memory weight parameters (e.g., local seed banks or persistent vegetative state).
#' @param Alfa An optional square numeric matrix of dimensions S x S defining interspecific
#'   competition coefficients between all pairs of species. If \code{NULL} (default), no interspecific competition is applied.
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
masterEqMetacomm <- function(Meta.pool,
                             d.spp,
                             FF = NULL,
                             m.pool,
                             Js,
                             M.migra,
                             id.fixed = NULL,
                             comm.fixed = NULL,
                             init.comm = NULL,
                             lottery = TRUE,
                             it = 100,
                             prop.dead.by.it = 0.05,
                             Ea = 1e-5,
                             Ts = 293.15,
                             m.temp = 0,
                             Alfa = NULL,
                             verbose = TRUE) {


  #### 1. VALIDATE INPUTS ####

  validateMetaInputs(Meta.pool = Meta.pool,
                     d.spp = d.spp,
                     FF = FF,
                     m.pool = m.pool,
                     Js = Js,
                     M.migra = M.migra,
                     id.fixed = id.fixed,
                     comm.fixed = comm.fixed,
                     init.comm = init.comm,
                     lottery = lottery,
                     it = it,
                     prop.dead.by.it = prop.dead.by.it,
                     Ea = Ea,
                     Ts = Ts,
                     m.temp = m.temp,
                     Alfa = Alfa
  )

  #### 2. INITIALIZATION AND DATA NORMALIZATION ####

  # Structural Dimension References
  S <- length(Meta.pool)
  C <- length(Js)

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


  #### 3. THERMAL DEPENDENCE CALCULATIONS ####

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


  #### 4. COALESCENT ASSEMBLY PHASE ####

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

    # Calculate interspecific competition overlap (this is skipped if Alfa is NULL)
    if (!is.null(Alfa)) {
      overlap <- Alfa %*% Meta
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

  #### 5. LOTTERY DYNAMICS ####
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
    generations <- seq(1, it, 1 / prop.dead.by.it)

    # Loop over iterations of the lottery dynamic
    for (iteration in 1:it) {

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

      # Re-calculate interspecific competition overlap (this is skipped if Alfa is NULL)
      if (!is.null(Alfa)) {
        overlap <- Alfa %*% Meta
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

      if(verbose){ cat("lottery iteration", iteration, "of", it, "\n")}
    }
  }

  return(Meta)
}
