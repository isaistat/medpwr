#' @title Sample Size for Log-Transformed t-Test (O’Keeffe et al., 2017)
#'
#' @description
#' Computes the required sample size for comparing two independent positively skewed outcome groups
#' using a log-transformed t-test, as described by O’Keeffe et al. (2017).
#' The method accounts for unequal standard deviations or ranges
#' and supports both one- and two-sided alternatives.
#'
#' @param m Numeric vector of length 2, representing the two group medians (or means) to be compared.
#' @param s Numeric vector of length 2, representing variability (standard deviation or range) for each group.
#' @param variability Character string specifying the type of variability measure:
#'   \itemize{
#'     \item `"SD"` (default) — values in `s` are standard deviations.
#'     \item `"range"` — values in `s` are ranges, which will be internally divided by 4 to approximate SDs.
#'   }
#' @param alternative Character string specifying the type of hypothesis test:
#'   \itemize{
#'     \item `"two.sided"` (default) — two-sided test.
#'     \item `"one.sided"` — one-sided test.
#'   }
#' @param alpha Numeric value for the significance level (default = 0.05).
#' @param power Numeric value for desired power (default = 0.80).
#'
#' @details
#' The function uses the log transformation to handle skewed data (e.g., median-based analysis)
#' and calculates the log-transformed standard deviations as:
#' \deqn{\sigma^2 = \log((1/2) + \sqrt{(1/4) + (s^2 / m^2)})}
#' The resulting formula for sample size per group is:
#' \deqn{n = \frac{(\sigma_1^2 + \sigma_2^2)(z_{\alpha} + z_{\beta})^2}{(\log(m_1) - \log(m_2))^2}}
#'
#' @return A list (invisible) containing:
#' \itemize{
#'   \item `n_per_group` — required sample size per group.
#'   \item `n_total` — total required sample size.
#'   \item `sigma` — log-transformed SDs for both groups.
#'   \item `log_diff` — log difference of the medians (effect size).
#' }
#'
#' @examples
#' # Example with SDs
#' med.test(m = c(30, 20), s = c(10, 10), variability = "SD", power = 0.8)
#'
#' # Example with ranges
#' med.test(m = c(25, 18), s = c(40, 35), variability = "range", alpha = 0.05)
#'
#' @references
#' O’Keeffe AG, Ambler G, Barber JA (2017). 
#' Sample size calculations based on a difference in medians for positively skewed outcomes in health care studies. 
#' BMC Med Res Methodol. 2017 Dec 2;17(1):157. doi: 10.1186/s12874-017-0426-1. PMID: 29197347; PMCID: PMC5712177.
#'
#' @export

med.test <- function(m = c(m1, m2),
                     s = c(s1, s2),
                     variability="SD",
                     alternative="two.sided",
                     alpha = 0.05,
                     power = 0.80)
{

  if (length(m) != 2 || length(s) != 2) {
    message("Error: Both m and phi must be vectors of length 2 (for two samples).")
    return(invisible("Error: Both m and phi must be vectors of length 2 (for two samples)."))
  }

  if (!variability %in% c("SD", "range")) {
    message("Error: variability must be either 'SD' or 'range'.")
    return(invisible("Error: variability must be either 'SD' or 'range'."))
  }
  if (!alternative %in% c("two.sided", "one.sided")) {
    message("Error: alternative must be either 'two.sided' or 'one.sided'.")
    return(invisible("Error: alternative must be either 'two.sided' or 'one.sided'."))
  }

  # --- Difference of logs ---
  t <- log(m[1]) - log(m[2])

  # Inside your function, after checking phi lengths
  if(s[1] == s[2]) {
    warning("s1=s2, provided. \n Using unequal variances can provide more precise sample size calculation.")
  }

  # --- Sigma calculation ---
  for (i in 1:2) {
    if (variability == "SD")
    {
      s_1=s[1]
      s_2=s[2]
    }
    else if (variability == "range")
    {
      s_1=s[1]/4
      s_2=s[2]/4
    }
  }
  # --- Log-transformed SDs as per O’Keeffe et al. ---
  sigma1 <- log(0.5 + sqrt(0.25 + (s_1^2 / m[1]^2)))
  sigma2 <- log(0.5 + sqrt(0.25 + (s_2^2 / m[2]^2)))

  # --- z-values ---
  z_alpha <- if(alternative!="two.sided") qnorm(1-alpha) else qnorm(1-(alpha/2))
  z_beta  <- qnorm(power)

  # --- Sample size per group ---
  n_per_group <- ((sigma1 + sigma2) * (z_alpha + z_beta)^2) / (t^2)
  n_total <- 2 * ceiling(n_per_group)

  # --- Output ---
  cat("\n---- Log-transformed t-test sample size (O’Keeffe 2017) ----\n")
  cat("*Total sample size =", n_total, "*\n*Each group requires =", ceiling(n_per_group), "*\n\n")
  cat("Effect size (log difference of medians) =", round(t, 4), "\n")
  cat("Level of significance =", alpha, "; Power =", round(power*100,0), "%\n")
  cat("Log-transformed SDs =", round(c(sigma1, sigma2), 4), "\n")
  cat("-------------------------------------------\n")

  invisible(list(
    n_per_group = ceiling(n_per_group),
    n_total = n_total,
    sigma = c(sigma1, sigma2),
    log_diff = t
  ))
}
