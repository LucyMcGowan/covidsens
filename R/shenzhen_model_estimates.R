#'  Model estimates from Shenzhen Bayesian Logistic Regression
#'
#' A list containing components of the Bayesian logistic model fit to the Shenzhen Data
#'
#' @format A list with 8 compoments, each with 5,000 observations for the 5,000 runs:
#' \describe{
#'   \item{beta_0}{The intercept from the logistic regression model}
#'   \item{beta_1}{Coefficient for the first component of the polynomial term in the logistic regression model}
#'   \item{beta_2}{Coefficient for the second component of the polynomial term in the logistic regression model}
#'   \item{beta_3}{Coefficient for the third component of the polynomial term in the logistic regression model}
#'   \item{mu}{Predicted values from the logistic regression model}
#'   \item{sens}{Sensitivity estimate from the logistic regression model}
#'   \item{log_lik}{Log likelihood from the logistic regression model}
#'   \item{lp__}{Log density up to a constant}
#' }
#' @source Zhang, Z. et al. Insights into the practical effectiveness of RT-PCR testing for SARS-CoV-2 from serologic data, a cohort study. doi:10.1101/2020.09.01.20182469.
#'
#' Kucirka, L. M., Lauer, S. A. & Laeyendecker, O. Variation in false-negative rate of reverse transcriptase polymerase chain reaction–based SARS-CoV-2 tests by time since exposure. Ann. Intern. Med. (2020)
"shenzhen_model_estimates"
