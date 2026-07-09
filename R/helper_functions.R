#' Random multinomial selection of new individuals based on mixed migration/regional probabilities
#'
#' @param n
#' @param dead.by.it
#' @param M.pool
#' @param m.pool
#'
#' @returns
#' @export
#'
#' @examples
born <- function(n, dead.by.it, M.pool, m.pool) {
  rmultinom(1, dead.by.it, (1 - m.pool) * (n / sum(n)) + m.pool * M.pool)
}

#' Handles specific subtraction adjustments during the mortality sub-phase
#'
#' @param n
#' @param change
#'
#' @returns
#' @export
#'
#' @examples
change <- function(n, change) {
  rmultinom(1, change, n)
}

#' Standardizes vector proportions alongside the global pool (kept for structural preservation)
#'
#' @param m
#' @param m.pool
#'
#' @returns
#' @export
#'
#' @examples
m_to_1 <- function(m, m.pool) {
  (1 - m.pool) * m / sum(m)
}
