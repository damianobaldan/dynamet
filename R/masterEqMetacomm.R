#' Simulate Metacommunity Dynamics via Coalescent Assembly and Lottery Phases
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
#' @param FF A numeric matrix of dimensions S x C representing local environmental
#'   filtering filters. Individual coefficients must range between 0 and 1, where 1 indicates
#'   perfect environmental match (no filtering penalty).
#' @param m.pool A single numeric value between 0 and 1 defining the probability of
#'   recruitment originating from the global regional pool rather than local/neighboring sources.
#' @param Js A numeric vector of length C setting the local carrying capacity (total
#'   individual slots) for each community patch.
#' @param M.migra A square numeric matrix of dimensions C x C establishing spatial
#'   migration connectivity and dispersal probabilities between patches. Cannot be NULL.
#' @param id.fixed An optional numeric vector containing indices of communities whose
#'   compositions are static and locked during simulation steps.
#' @param comm.fixed An optional numeric vector of length S outlining the fixed relative
#'   abundance profile assigned to the patches listed in \code{id.fixed}.
#' @param Lottery A single logical value. If \code{TRUE}, triggers the demographic turnover
#'   iterative lottery phase after initial community assembly.
#' @param it A single numeric integer specifying the total timeline steps/iterations to
#'   run inside the lottery loop.
#' @param prop.dead.by.it A single numeric fraction (0, 1) defining the baseline mortality
#'   turnover rate, anchored at the coldest patch.
#' @param id.obs An optional numeric vector of community indices reserved for observation,
#'   validation, or structural tracking.
#' @param Ea A single numeric value representing the activation energy (in eV) utilized to
#'   scale temperature-driven mortality.
#' @param Ts A numeric vector of length C specifying the local patch ambient temperatures
#'   expressed strictly in Kelvin units (> 0).
#' @param m.temp An optional numeric vector of length C, or a matrix of dimensions S x C,
#'   dictating community memory weight parameters (e.g., local seed banks or persistent vegetative state).
#' @param Alfa A square numeric matrix of dimensions S x S defining interspecific
#'   competition coefficients between all pairs of species.
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
#' @examples
#' \dontrun{
#' # Assuming inputs S = 3 species, C = 2 communities:
#' final_meta <- masterEqMetacomm(
#'   Meta.pool = c(0.5, 0.3, 0.2), d.spp = c(1, 1, 1),
#'   FF = matrix(1, nrow = 3, ncol = 2), m.pool = 0.1, Js = c(100, 150),
#'   M.migra = diag(2), Lottery = TRUE, it = 10, prop.dead.by.it = 0.1,
#'   Ea = 0.65, Ts = c(293.15, 298.15), m.temp = NULL, Alfa = diag(3)
#' )
#' }
masterEqMetacomm <- function(Meta.pool,
                             d.spp,
                             FF,
                             m.pool,
                             Js,
                             M.migra = NULL,
                             id.fixed = NULL,
                             comm.fixed = NULL,
                             Lottery,
                             it,
                             prop.dead.by.it,
                             id.obs,
                             Ea,
                             Ts,
                             m.temp,
                             Alfa,
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
                     Lottery = Lottery,
                     it = it,
                     prop.dead.by.it = prop.dead.by.it,
                     id.obs = id.obs,
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
    comm.fixed <- comm.fixed / sum(comm.fixed)
  }

  # Initialize community memory tracking (m.temp)
  if (is.null(m.temp)) {
    m.temp <- rep(0, C)
  } else if (!is.matrix(m.temp) && length(m.temp) == 1) {
    m.temp <- rep(m.temp, C)
  }



  #### 3. THERMAL DEPENDENCE CALCULATIONS ####

  # Estimate standard constant based on the minimum temperature patch to anchor baseline mortality
  # 8.62e-5 represents the Boltzmann constant in eV/K
  min.dead.Tmin <- prop.dead.by.it / exp(-Ea / (min(Ts) * 8.62e-5))

  # Scale individual community mortality rates based on local temperature (Arrhenius relationship)
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

  #### 4. COALESCENT ARGUMENT ####

  # Seed each community with exactly 1 individual based on regional pool and filters
  Meta <- matrix(NA, nrow = length(Meta.pool), ncol = ncol(M.migra))
  for (i in 1:ncol(M.migra)) {
    Meta[,i] <- rmultinom(1, 1, Meta.pool * d.spp * FF[,i])
  }

  # Loop over communities to perform the coalescent assembly (until max(Js) is reached) Iteratively scale up communities until they reach their respective carrying capacity (Js)
  for (ii in 2:max(Js)) {

    # Track communities that have not yet reached their carrying capacity
    id.j <- which(Js >= ii)

    # Keep fixed communities scaled to current global abundance level
    Meta[, id.fixed] <- comm.fixed * (ii - 1)

    # Calculate potential recruits based on local abundance and spatial migration (apply dispersal and env. filter)
    Pool.neighbor <- (Meta %*% M.migra)  * d.spp  * FF
    # Pool.neighbor <- Pool.neighbor * d.spp          # Apply species dispersal constraints
    # Pool.neighbor <- Pool.neighbor * FF             # Apply local environmental filters

    # Calculate interspecific competition overlap and scale by column
    overlap       <- Alfa %*% Meta
    overlap       <- sweep(overlap, 2, colSums(overlap), FUN = "/")
    Pool.neighbor <- Pool.neighbor * (1 - overlap)  # Reduce recruitment pool by competitive pressure

    # Assign new individuals to available spaces
    if (length(id.j) > 1) {
      new <- apply(Pool.neighbor[, id.j], 2, born, dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool)
      Meta[, id.j] <- Meta[, id.j] + new
    } else {
      Meta[, id.j] <- Meta[, id.j] + born(n = Pool.neighbor[, id.j], dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool)
    }


    # Print iteration
    if(verbose){ cat("Coalescent construction in J:", ii, "of", max(Js), "\n") }

  }

  #### 5. LOTTERY DYNAMICS ####
  if (Lottery == TRUE) {

    # Mirror state for community memory tracking
    Meta.lag         <- Meta

    # Set fixed communities to maximum carrying capacity
    Meta[, id.fixed] <- round(comm.fixed * max(Js), 0)

    # Generate timeline checkpoints
    generations      <- seq(1, it, 1 / prop.dead.by.it)

    # Loop over iterations of the lottery dynamic
    for (iteration in 1:it) {

      # update community memory cache if at checkpoint
      if (any(generations==it)) { Meta.lag <- Meta }

      # --- Death Sub-phase ---
      for (dead in 1:max.dead.by.it) {

        # Target patches experiencing active mortality
        id.dead <- which(dead.by.it >= dead)

        # Select and remove individuals (poorly filtered species have higher mortality risks)
        if (length(id.dead) > 1) {
          Meta[, id.dead] <- Meta[, id.dead] - apply(Meta[, id.dead] * (1.001 - FF[, id.dead]), 2, FUN = change, change = 1)
        }
        if (length(id.dead) == 1) {
          Meta[, id.dead] <- Meta[, id.dead] - change(n = Meta[, id.dead] * (1.001 - FF[, id.dead]), change = 1)
        }
      }

      # --- Recruitment Sub-phase ---

      # Calculate potential recruits based on local abundance and spatial migration (apply dispersal and env. filter)
      Pool.neighbor <- (Meta %*% M.migra)  * d.spp  * FF
      # Pool.neighbor <- Pool.neighbor * d.spp          # Apply species dispersal constraints
      # Pool.neighbor <- Pool.neighbor * FF             # Apply local environmental filters

      # Factor in community memory ("seed bank" preservation vs new migrants)
      Pool.neighbor <- Pool.neighbor * (1 - m.temp) + Meta.lag * (m.temp)

      # Re-evaluate competitive structural overlap
      overlap       <- Alfa %*% Meta
      overlap       <- sweep(overlap, 2, colSums(overlap), FUN = "/")
      Pool.neighbor <- Pool.neighbor * (1 - overlap)

      # Normalize Pool.neighbor by column and add the contribution of m.pool
      Pool.norm <- sweep(Pool.neighbor, 2, colSums(Pool.neighbor), FUN = "/")
      Prob.mat  <- (1 - m.pool) * Pool.norm + m.pool * Meta.pool

      # Fill vacancies using multinomial sampling across local, neighboring, and regional pools
      for (i in seq_along(dead.by.it)) {
        Meta[, i] <- Meta[, i] + rmultinom(1, size = dead.by.it[i], prob = Prob.mat[, i])
      }

      # print iteration
      if(verbose){ cat("Lottery iteration", iteration, "of", it, "\n")}

    }
  }

  return(Meta)

}
