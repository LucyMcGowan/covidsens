#' Get test sensitivity x days from crossing the threshold of detection
#'
#' @param time Numeric. Number of days from crossing the threshold of detection
#' @param shape Numeric. Shape of the gamma distribution used to estimate
#'   test sensitivity. Default is 2.32
#' @param rate Numeric. Rate of the gamma distribution used to estimate test
#'   sensitivity. Default is 0.23
#'
#' @export
#'
#' @examples
#' get_sensitivity(1)
get_sensitivity <- function(time, shape = 2.32, rate = 0.23) {
  stats::dgamma(time, shape, rate) /
    stats::dgamma((shape - 1) / rate, shape, rate) * 0.99
}

#' Get test false negative rate x days from crossing the threshold of detection
#'
#' @param time Numeric. Number of days from crossing the threshold of detection
#' @param shape Numeric. Shape of the gamma distribution used to estimate
#'   test sensitivity. Default is 2.33
#' @param rate Numeric. Rate of the gamma distribution used to estimate test
#'   sensitivity. Default is 0.24
#'
#' @export
#'
#' @examples
#' get_fnr(1)
get_fnr <- function(time, shape = 2.33, rate = 0.24) {
  1 - get_sensitivity(time, shape, rate)
}
