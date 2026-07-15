#' Random Multinomial Selection of New Recruits
#'
#' @description
#' Simulates the recruitment of new individuals into community vacancies using a
#' multinomial draw. The sampling probabilities are calculated as a weighted blend
#' of local/neighboring species abundances and the regional species pool.
#'
#' @param n A numeric vector representing the raw species recruitment weights or
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
born <- function(n, dead.by.it, M.pool, m.pool) {
  stats::rmultinom(1, dead.by.it, (1 - m.pool) * (n / sum(n)) + m.pool * M.pool)
}



#' Handle Abundance Reductions During Mortality Sub-Phases
#'
#' @description
#' Executes a multinomial draw to determine which specific individuals die within a
#' community patch. The probability of death is weighted by current species
#' abundances and their respective environmental filtering penalties.
#'
#' @param n A numeric vector representing the weighted vulnerability or abundance
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
change <- function(n, change) {
  stats::rmultinom(1, change, n)
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

