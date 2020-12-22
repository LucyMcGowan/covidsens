#' Extract a single estimate of the Shenzhen sensitivity based on the Bayesian
#' logistic regression
#'
#' @param x Numeric. Days since symptom onset
#'
#' @export
#'
#' @examples
#' get_shenzhen_estimate(0)
get_shenzhen_estimate <- function(x) {

  ## recreate polynomial object
  pol <- list()
  attr(pol, "coefs") <- list(alpha = c(13.71429, 12.59599, 13.53112),
                             norm2 = c(1, 21, 1020.286, 55471.585, 2413933.056))
  attr(pol, "degree") <- c(1, 2, 3)
  class(pol) <- c("poly")

  id <- sample(1:5000, 1)
  x <- x + 9
  if (x == Inf) {
    return(1)
  }
  if (x < - 12) {
    return(1)
  }
  p <- stats::predict(pol, x)
  1 - inv_logit(covidsens::shenzhen_model_estimates$beta_0[id] +
                  covidsens::shenzhen_model_estimates$beta_1[id] * p[,1] +
                  covidsens::shenzhen_model_estimates$beta_2[id] * p[,2] +
                  covidsens::shenzhen_model_estimates$beta_3[id] * p[,3])
}

inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}



