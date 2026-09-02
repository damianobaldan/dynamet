#' Validate Inputs for Mutualistic Metacommunity Models
#'
#' Checks the dimensions, types, and constraints of input arguments before
#' running Master Equation or simulation models for bipartite mutualistic networks.
#'
#' @inheritParams masterEqMutualisticMetacomm
#'
#' @return Invisibly returns `TRUE` if all inputs are valid. Throws an error otherwise.
#' @export
validateMutualisticInputs <- function(Meta.pool, m.pool, Js.plant, Js.animal,
                                      filter.env, dispersal.spp, Q.LH,
                                      M.dist, D50, m.max,
                                      id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                                      it, niches.mean, niches.sd, Sp, Sa) {

  S <- Sp + Sa
  N_comm <- nrow(M.dist)

  # Structural length & dimension checks
  if (length(Meta.pool) != S) {
    stop("Length of 'Meta.pool' must equal Sp + Sa (", S, ").")
  }
  if (length(Js.plant) != N_comm) {
    stop("Length of 'Js.plant' (", length(Js.plant), ") must match number of communities in 'M.dist' (", N_comm, ").")
  }
  if (length(Js.animal) != N_comm) {
    stop("Length of 'Js.animal' (", length(Js.animal), ") must match number of communities in 'M.dist' (", N_comm, ").")
  }
  if (nrow(filter.env) != S || ncol(filter.env) != N_comm) {
    stop("Dimensions of 'filter.env' (", nrow(filter.env), "x", ncol(filter.env),
         ") must be (Sp + Sa) x N_comm (", S, "x", N_comm, ").")
  }
  if (length(dispersal.spp) != S) {
    stop("Length of 'dispersal.spp' must equal Sp + Sa (", S, ").")
  }
  if (length(niches.mean) != S || length(niches.sd) != S) {
    stop("Lengths of 'niches.mean' and 'niches.sd' must equal Sp + Sa (", S, ").")
  }
  if (nrow(M.dist) != ncol(M.dist)) {
    stop("'M.dist' must be a square distance matrix.")
  }

  # Fixed community parameter checks
  if (length(id.fixed) > 0) {
    if (any(id.fixed < 1 | id.fixed > N_comm)) {
      stop("All indices in 'id.fixed' must be between 1 and the number of communities (", N_comm, ").")
    }
    if (is.null(comm.fixed) || length(comm.fixed) != S) {
      stop("'comm.fixed' must be provided with length equal to Sp + Sa (", S, ") when 'id.fixed' is non-empty.")
    }
    if (is.null(D50.fixed) || is.null(m.max.fixed)) {
      stop("'D50.fixed' and 'm.max.fixed' must be specified when 'id.fixed' is non-empty.")
    }
  }

  # Range and domain constraints
  if (m.pool < 0 || m.pool > 1) {
    stop("'m.pool' must be a probability between 0 and 1.")
  }
  if (any(filter.env < 0)) {
    stop("All values in 'filter.env' must be non-negative.")
  }
  if (any(apply(filter.env, 2, sum) == 0)) {
    stop("At least one community in 'filter.env' has all species performance set to zero.")
  }
  if (any(niches.sd <= 0)) {
    stop("All values in 'niches.sd' must be strictly positive (> 0).")
  }
  if (it <= 0 || it != as.integer(it)) {
    stop("'it' must be a positive integer.")
  }

  return(invisible(TRUE))
}


#' Stochastic Master Equation Engine for Bipartite Mutualistic Metacommunities
#'
#' Simulates the stochastic abundance dynamics of a plant-animal metacommunity over time
#' using a Master Equation framework based on spatial dispersal, environmental filtering,
#' density-dependent competition, and mutualistic interaction networks.
#'
#' @param Meta.pool Numeric vector. Regional species pool abundances (Plants + Animals).
#' @param m.pool Numeric. Global immigration rate from the regional pool (0 to 1).
#' @param Js.plant Numeric vector. Local community carrying capacities for plants.
#' @param Js.animal Numeric vector. Local community carrying capacities for animals.
#' @param filter.env Matrix. Abiotic environmental filter values (Species x Communities).
#' @param dispersal.spp Numeric vector. Relative dispersal abilities of all species.
#' @param Q.LH Numeric. Scaling exponent for mutualistic benefit.
#' @param M.dist Matrix. Symmetric spatial distance matrix between local communities.
#' @param D50 Numeric. Half-distance dispersal decay parameter.
#' @param m.max Numeric. Maximum local migration rate between neighboring sites.
#' @param id.fixed Integer vector. Community indices with fixed abundances. Default is `c()`.
#' @param D50.fixed Numeric. Half-distance decay parameter for fixed communities. Default is `NULL`.
#' @param m.max.fixed Numeric. Maximum migration rate for fixed communities. Default is `NULL`.
#' @param comm.fixed Numeric vector. Relative species composition for fixed communities. Default is `NULL`.
#' @param it Integer. Total number of iterations to run. Default is `100`.
#' @param niches.mean Numeric vector. Trait means for all species.
#' @param niches.sd Numeric vector. Trait standard deviations for all species.
#' @param Sp Integer. Total number of plant species.
#' @param Sa Integer. Total number of animal species.
#'
#' @details
#' The model evaluates interaction matrices deterministically at each timestep based
#' on niche overlap integration between species, taking into account current abundances.
#' New species states for each local community are generated via multinomial draws
#' driven by the combined effects of immigration and local carrying capacities (\code{Js.plant}
#' and \code{Js.animal}).
#'
#' @return An integer matrix of dimension \code{(Sp + Sa) x N_comm} representing final species
#'   abundances per community at the end of \code{it} iterations.
#' @export
#'
#' @seealso \code{\link{validateMutualisticInputs}}
masterEqMutualisticMetacomm <- function(Meta.pool, m.pool, Js.plant, Js.animal,
                                        filter.env, dispersal.spp, Q.LH,
                                        M.dist, D50, m.max,
                                        id.fixed = c(), D50.fixed = NULL,
                                        m.max.fixed = NULL, comm.fixed = NULL,
                                        it = 100, niches.mean, niches.sd,
                                        Sp, Sa) {

  # Validate inputs prior to execution
  validateMutualisticInputs(Meta.pool = Meta.pool, m.pool = m.pool,
                            Js.plant = Js.plant, Js.animal = Js.animal,
                            filter.env = filter.env, dispersal.spp = dispersal.spp,
                            Q.LH = Q.LH, M.dist = M.dist, D50 = D50, m.max = m.max,
                            id.fixed = id.fixed, D50.fixed = D50.fixed,
                            m.max.fixed = m.max.fixed, comm.fixed = comm.fixed,
                            it = it, niches.mean = niches.mean, niches.sd = niches.sd,
                            Sp = Sp, Sa = Sa)


  # helper functions
  migration.matrix.kernel.all <- function(M.dist, m.pool, D50, m.max, id.fixed, D50.fixed, m.max.fixed) {
    diag(M.dist) <- 0
    M.dist <- M.dist - min(M.dist, na.rm = TRUE)
    b <- -log(0.5) / D50
    M.migra <- m.max * exp(-b * M.dist)

    if (length(id.fixed) > 0) {
      b.fixed <- -log(0.5) / D50.fixed
      M.migra[id.fixed, ] <- m.max.fixed * exp(-b.fixed * M.dist[id.fixed, ])
    }
    diag(M.migra) <- 1
    M.migra <- apply(M.migra, 2, m_to_1, m.pool)
    return(M.migra)
  }

  m_to_1 <- function(m, m.pool) (1 - m.pool) * m / sum(m)

  S <- Sp + Sa
  N_comm <- nrow(M.dist)

  # Normalize global regional pool relative abundances
  Meta.pool.plant <- Meta.pool[1:Sp] / sum(Meta.pool[1:Sp])
  Meta.pool.animal <- Meta.pool[(Sp + 1):S] / sum(Meta.pool[(Sp + 1):S])
  Meta.pool_norm <- c(Meta.pool.plant, Meta.pool.animal)

  # Distance-based dispersal/migration kernel
  M.migra <- migration.matrix.kernel.all(M.dist = M.dist, m.pool = m.pool,
                                         D50 = D50, m.max = m.max,
                                         id.fixed = id.fixed, D50.fixed = D50.fixed,
                                         m.max.fixed = m.max.fixed)

  # Format fixed community proportions if provided
  if (length(id.fixed) > 0 && !is.null(comm.fixed)) {
    comm.fixed[1:Sp] <- comm.fixed[1:Sp] / sum(comm.fixed[1:Sp])
    comm.fixed[(Sp + 1):S] <- comm.fixed[(Sp + 1):S] / sum(comm.fixed[(Sp + 1):S])
  }

  # Initialize metacommunity state matrix via multinomial draws from regional pool
  Meta <- matrix(0, nrow = S, ncol = N_comm)
  for (j in 1:N_comm) {
    if (j %in% id.fixed) {
      Meta[1:Sp, j] <- round(comm.fixed[1:Sp] * Js.plant[j])
      Meta[(Sp + 1):S, j] <- round(comm.fixed[(Sp + 1):S] * Js.animal[j])
    } else {
      Meta[1:Sp, j] <- rmultinom(1, Js.plant[j], Meta.pool.plant)
      Meta[(Sp + 1):S, j] <- rmultinom(1, Js.animal[j], Meta.pool.animal)
    }
  }

  # Static trait matrices
  Traits.mu <- matrix(niches.mean, byrow = FALSE, ncol = N_comm, nrow = S)
  Traits.sd <- matrix(niches.sd, byrow = FALSE, ncol = N_comm, nrow = S)

  mutualist.overlap <- matrix(0, nrow = S, ncol = N_comm)
  competition.overlap <- matrix(0, nrow = S, ncol = N_comm)

  # Fixed Stochastic Master Loop
  for (step in 1:it) {

    # Calculate density-dependent interaction overlaps based on current state
    for (local in 1:N_comm) {
      ov.all <- net.from.N.ov(Ns = Meta[, local], niches.mean = Traits.mu[, local],
                              niches.sd = Traits.sd[, local], Sp = Sp, Sa = Sa)

      # Mutualism overlap
      vector.ov.mutual.plant <- apply(ov.all[[1]] + 1, 1, sum) / sum(ov.all[[1]] + 1) + 1e-7
      vector.ov.mutual.animal <- apply(ov.all[[1]] + 1, 2, sum) / sum(ov.all[[1]] + 1) + 1e-7
      vector.ov.mutual.plant <- (vector.ov.mutual.plant / max(vector.ov.mutual.plant)) / Traits.sd[1:Sp, local]^2
      vector.ov.mutual.animal <- (vector.ov.mutual.animal / max(vector.ov.mutual.animal)) / Traits.sd[(1 + Sp):S, local]^2

      mutualist.overlap[, local] <- c(vector.ov.mutual.plant, vector.ov.mutual.animal)

      # Competition overlap
      vector.ov.commp.plant <- apply(1 - ov.all[[2]], 1, sum, na.rm = TRUE) / sum(1 - ov.all[[2]], na.rm = TRUE) + 1e-7
      vector.ov.commp.plant <- vector.ov.commp.plant / max(vector.ov.commp.plant)
      vector.ov.commp.animal <- apply(1 - ov.all[[3]], 1, sum, na.rm = TRUE) / sum(1 - ov.all[[3]], na.rm = TRUE) + 1e-7
      vector.ov.commp.animal <- vector.ov.commp.animal / max(vector.ov.commp.animal)

      rel.competition <- c(vector.ov.commp.plant, vector.ov.commp.animal)
      competition.overlap[, local] <- pmax(0, rel.competition)
    }

    # Effective immigration pool weighting dispersal, mutualism, competition, and abiotic filtering
    Pool.neighbor <- (Meta %*% M.migra) * (mutualist.overlap^Q.LH) * competition.overlap * filter.env * dispersal.spp

    # Stochastic Multinomial Sampling across non-fixed local communities
    for (j in 1:N_comm) {
      if (j %in% id.fixed) {
        Meta[1:Sp, j] <- round(comm.fixed[1:Sp] * Js.plant[j])
        Meta[(Sp + 1):S, j] <- round(comm.fixed[(Sp + 1):S] * Js.animal[j])
      } else {
        # Plant Multinomial Transition
        plant_pool <- Pool.neighbor[1:Sp, j]
        prob_plant <- (1 - m.pool) * (plant_pool / sum(plant_pool)) + m.pool * Meta.pool_norm[1:Sp]
        Meta[1:Sp, j] <- rmultinom(1, size = Js.plant[j], prob = prob_plant)

        # Animal Multinomial Transition
        animal_pool <- Pool.neighbor[(Sp + 1):S, j]
        prob_animal <- (1 - m.pool) * (animal_pool / sum(animal_pool)) + m.pool * Meta.pool_norm[(Sp + 1):S]
        Meta[(Sp + 1):S, j] <- rmultinom(1, size = Js.animal[j], prob = prob_animal)
      }
    }
  }

  rownames(Meta) <- c(paste0("Plant_", 1:Sp), paste0("Animal_", 1:Sa))
  colnames(Meta) <- paste0("Comm_", 1:N_comm)

  return(Meta)
}
