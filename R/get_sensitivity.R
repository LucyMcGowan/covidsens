#' Get test sensitivity x days from crossing the threshold of detection
#'
#' @param time Numeric. Number of days from crossing the threshold of detection
#' @param shape Numeric. Shape of the gamma distribution used to estimate
#'   test sensitivity. Default is 2.33
#' @param rate Numeric. Rate of the gamma distribution used to estimate test
#'   sensitivity. Default is 0.24
#' @param max_sensitivity Numeric. Value between 0 and 1 indicating the maximum
#'   sensitivity. Default is 0.99.
#'
#' @export
#'
#' @examples
#' get_sensitivity(1)
get_sensitivity <- function(time, shape = 2.33, rate = 0.24, max_sensitivity = 0.99) {
  stats::dgamma(time, shape, rate) /
    stats::dgamma((shape - 1) / rate, shape, rate) * max_sensitivity
}

#' Get test false negative rate x days from crossing the threshold of detection
#'
#' @param time Numeric. Number of days from crossing the threshold of detection
#' @param shape Numeric. Shape of the gamma distribution used to estimate
#'   test sensitivity. Default is 2.33
#' @param rate Numeric. Rate of the gamma distribution used to estimate test
#'   sensitivity. Default is 0.24
#' @param max_sensitivity Numeric. Value between 0 and 1 indicating the maximum
#'   sensitivity. Default is 0.99.
#'
#' @export
#'
#' @examples
#' get_fnr(1)
get_fnr <- function(time, shape = 2.33, rate = 0.24, max_sensitivity = 0.99) {
  1 - get_sensitivity(time, shape, rate, max_sensitivity)
}
