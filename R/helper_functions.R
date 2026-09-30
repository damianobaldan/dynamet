#' Random Multinomial Selection of New Recruits
#'
#' @description
#' Simulates the recruitment of new individuals into community vacancies using a
#' multinomial draw. The sampling probabilities are calculated as a weighted blend
#' of local/neighboring species abundances and the regional species pool.
#'
#' @param probs A numeric vector representing the raw species recruitment weights or
#'   abundances derived from local and neighboring patches.
#' @param dead.by.it An integer specifying the number of individuals (vacancies)
#'   to sample for recruitment during this step.
#' @param M.pool A numeric vector of species relative abundances within the global
#'   regional pool (\code{Meta.pool}).
#' @param m.pool A single numeric value between 0 and 1 representing the probability
#'   that a recruit originates from the regional pool rather than local sources.
#'
#' @return A numeric matrix or vector containing the counts of recruited individuals
#'   assigned to each species.
#'
#' @importFrom stats rmultinom
#' @keywords internal
#'
born <- function(probs, dead.by.it, M.pool, m.pool) {
  stats::rmultinom(1, dead.by.it, (1 - m.pool) * (probs / sum(probs)) + m.pool * M.pool)
}



#' Handle Abundance Reductions During Mortality Sub-Phases
#'
#' @description
#' Executes a multinomial draw to determine which specific individuals die within a
#' community patch. The probability of death is weighted by current species
#' abundances and their respective environmental filtering penalties.
#'
#' @param probs A numeric vector representing the weighted vulnerability or abundance
#'   profile of species within the targeted patch.
#' @param change An integer specifying the number of individuals to select and
#'   remove (kill) from the patch.
#'
#' @return A numeric matrix or vector containing the counts of individuals to be
#'   subtracted per species.
#'
#' @importFrom stats rmultinom
#' @keywords internal
#'
change <- function(probs, change) {
  stats::rmultinom(1, change, probs)
}

#' Prune Local Communities to Target Carrying Capacities
#'
#' Filters and downsamples a community abundance matrix when local site abundances
#' exceed specified target carrying capacities (e.g., during environmental disturbance).
#' Downsampling is performed via multinomial sampling, either neutrally based on current
#' abundances or weighted by local species fitness values.
#'
#' @param Meta A numeric matrix of species abundances where rows represent
#'   species (\eqn{S}) and columns represent local sites/cells (\eqn{C}).
#' @param Js A numeric vector of length \eqn{C} specifying target carrying
#'   capacities for each site.
#' @param FF An optional numeric matrix of species fitness or environmental
#'   suitability (\eqn{S \times C}). If provided, sampling probabilities are weighted
#'   by \code{Meta * FF}. Defaults to \code{NULL} (neutral culling).
#'
#' @details
#' For each site \eqn{c}, if total current abundance \eqn{\sum_i N_{i,c}} exceeds
#' \code{Js[c]}, individual organisms are downsampled to match the target capacity.
#' If \code{Js[c] <= 0}, all species abundances at site \eqn{c} are set to 0.
#'
#' @return A numeric matrix (\eqn{S \times C}) of filtered species abundances.
#'
#' @keywords internal
#'
filter_community <- function(Meta, Js, FF = NULL) {

  # Get parameters
  S <- nrow(Meta)
  C <- ncol(Meta)

  # Initialize filtered community matrix
  Meta_filtered <- Meta

  for (c in 1:C) {
    current_N <- sum(Meta[, c])
    target_J  <- Js[c]

    # Downsample only if current population exceeds new capacity
    if (current_N > target_J) {
      if (target_J <= 0) {
        Meta_filtered[, c] <- 0
      } else {
        # Calculate multinomial probabilities (fitness-weighted vs. neutral)
        if (!is.null(FF)) {
          weights <- Meta[, c] * FF[, c]
        } else {
          weights <- Meta[, c]
        }

        # Draw target_J surviving individuals from current occupants
        if (sum(weights) > 0) {
          Meta_filtered[, c] <- as.vector(stats::rmultinom(1, size = target_J, prob = weights ))
        } else {
          Meta_filtered[, c] <- 0
        }
      }
    }
  }

  return(Meta_filtered)
}



#' Standardize Local Abundance Proportions Relative to Global Pool
#'
#' @description
#' Downscales and normalizes local/neighboring recruitment weights to ensure that
#' the total available probability space accommodates the global immigration
#' contribution (\code{m.pool}).
#'
#' @param m A numeric vector of raw species relative abundances or local migration weights.
#' @param m.pool A single numeric value between 0 and 1 indicating the probability
#'   of recruitment originating from the regional pool.
#'
#' @return A numeric vector representing the adjusted, normalized probabilities of
#'   species recruitment originating strictly from the local neighborhood.
#'
#' @importFrom stats rmultinom
#' @keywords internal
#'
m_to_1 <- function(m, m.pool) {
  (1 - m.pool) * m / sum(m)
}



#' Get default parameters for metacommunity simulations
#'
#' @returns A named list of default simulation parameters.
#' @keywords internal
#'
getSimulationDefaults <- function() {
  list(
    d.spp = NULL,
    FF = NULL,
    alpha = NULL,
    init.comm = NULL,
    id.fixed = NULL,
    comm.fixed = NULL,
    prop.dead.by.it = 0.05,
    Ea = 1e-5,
    Ts = 293.15,
    m.temp = 0,
    lottery = TRUE,
    nIterations = 100,
    verbose = TRUE
  )
}

