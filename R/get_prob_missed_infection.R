#' Calculate the probability of a missed infection given a test + quarantine
#' strategy
#'
#' @param test_time Numeric. Date post-exposure test sample is collected
#' @param additional_quarantine_time Numeric. Additional days quarantined after
#'   test sample is collected. Default is 0
#' @param shape_fnr Numeric. Shape of the gamma distribution used to estimate
#'   test sensitivity. Default is 2.33
#' @param rate_fnr Numeric. Rate of the gamma distribution used to estimate test
#'   sensitivity. Default is 0.24
#' @param max_sensitivity Numeric. Value between 0 and 1 indicating the maximum
#'   sensitivity. Default is 0.99.
#' @param shape_exposure_to_threshold Numeric. Shape of the gamma distribution
#'   that describes the time from exposure to crossing the threshold of
#'   detection. Default is 1.98
#' @param shape_threshold_to_symptoms Numeric. Shape of the gamma distribution
#'   that describes the time from crossing the threshold of
#'   detection to symptom onset. Default is 2.16
#' @param rate Numeric. Rate of the gamma distributions that describe the time
#'   from exposure to crossing the threshold of detection and the time from
#'   crossing the threshold of detection to symptom onset. Default is 0.72.
#'
#' @export
#'
#' @examples
#' ## test sample collected 5 days post-exposure, quarantine until day 5
#' get_prob_missed_infection(5)
#'
#' ## test sample collected 7 days post-exposure, quarantine until day 9
#' get_prob_missed_infection(7, 2)
get_prob_missed_infection <- function(test_time,
                                      additional_quarantine_time = 0,
                                      shape_fnr = 2.33,
                                      rate_fnr = 0.24,
                                      max_sensitivity = 0.99,
                                      shape_exposure_to_threshold = 1.98,
                                      shape_threshold_to_symptoms = 2.16,
                                      rate = 0.72) {
  stats::integrate(top, 0, Inf, t = test_time, s = additional_quarantine_time,
                   shape_a = shape_fnr,
                   rate_a = rate_fnr,
                   max_sensitivity = max_sensitivity,
                   shape_y = shape_threshold_to_symptoms,
                   shape_x = shape_exposure_to_threshold,
                   rate = rate)$value
}

top <- function(x, t, s, shape_a, rate_a, max_sensitivity, shape_y, shape_x, rate) {
  (purrr::map_dbl(t - x, get_fnr, shape_a, rate_a, max_sensitivity) *
     (1 - G((t + s) - x, shape_y, rate)) * f(x, shape_x, rate))
}

G <- function(y, shape, rate) {
  stats::pgamma(y, shape, rate)
}

f <- function(x, shape, rate) {
  stats::dgamma(x, shape, rate)
}
