#' @title Estimated parameters for a Poisson INAR(1)
#' @name estimator_inar1
#'
#' @rdname estimator_inar1
#'
#' @usage yw(series)
#' sqd_diff(series)
#' kendall(series)
#'
#' @param series vector of (non-negative) integer time series
#'
#' @references Marcelo Bourguignon & Klaus L. P. Vasconcellos (2018).
#' The effects of additive outliers in INAR(1) process and robust estimation.
#'
#' @examples
#' # With the arima.sim() function
#' # this simulation can result in a negative lambda, since it's not a
#' # (non-negative) integer time series
#' yw(arima.sim(n = 1000, list(order = c(1,0,0), ar = 0.4)))
#' sqd_diff(arima.sim(n = 500, list(order = c(1,0,0), ar = 0.8)))
#' kendall(arima.sim(n = 300, list(order = c(1,0,0), ar = 0.9)))
#'
#' # With the simulation function
#' yw(sim_pois(1000, 0.4, 1))
#' sqd_diff(sim_pois(500, 0.4, 1))
#' kendall(sim_pois(300, 0.4, 1))
#' @export
yw <- function(series){ # order=1 (p=1)
  # Tamanho da serie tempora
  T_s <- length(series)
  # Media da serie temporal
  y_bar <- mean(series)

  rho <- {
    sum( (series[1:(T_s - 1)] - y_bar) * (series[2:T_s] - y_bar) ) /
      sum((series[1:T_s] - y_bar) ^ 2)
  }

  # Alpha estimado
  if(rho < 0){
    alpha <- 0
  } else{
    alpha <- rho
  }
  # Lambda estimado
  lambda <- (1 - alpha) * y_bar

  # Lista com as informacoes para a saida
  res <- list(call = match.call(),
              method = "yule-walker",
              order = 1,
              pacf = rho,
              coef = round(alpha, 4),
              lambda = lambda)

  return(res)
}

#' @rdname estimator_inar1
#' @usage NULL
#' @export
sqd_diff <- function(series){ # order=1 (p=1)
  # Tamanho da serie temporal
  T_s <- length(series)

  # Media da serie temporal
  y_bar <- mean(series)

  rho <- {
    sum( (series[1:(T_s - 1)] - y_bar) * (series[2:T_s] - y_bar) ) /
      sum((series[1:T_s] - y_bar) ^ 2)
  }

  # Lambda estimado
  lambda <- sum( (series[2:T_s] - series[1:(T_s - 1)])^2 ) / (2 * (T_s - 1))
  # Alpha estimado
  a <- 1 - (lambda / y_bar)
  if(a < 0){
    alpha <- 0
  } else{
    alpha <- a
  }

  # Lista com as informacoes para a saida
  res <- list(call = match.call(),
              method = "squared difference",
              order = 1,
              pacf = rho,
              coef = round(alpha, 4),
              lambda = lambda)

  return(res)
}

#' @rdname estimator_inar1
#' @usage NULL
#' @export
kendall <- function(series){ # order=1 (p=1)
  # Tamanho da serie temporal
  T_s <- length(series)

  # Media da serie temporal
  y_bar <- mean(series)

  # Estatistica K (Kendall, 1938)
  kendall_K <- function(series, len){
    K = 0
    # se nao subtrair 1 (order) mais uma vez, devolve NA quando chegar no
    # final da serie, pq ele ultrapassa o tamanho do vetor
    for(i in 1:(len - 1 - 1)){
      for(j in (i + 1):(len - 1)){
        q <- (series[j + 1] - series[i + 1]) * (series[j] - series[i])

        if(q > 0){
          K <- K + 1
        } else if(q == 0){
          K <- K
        } else{
          K <- K - 1
        }
      }
    }

    return(K)
  }
  K <- kendall_K(series, T_s)

  rho <- 2 * K / (T_s * (T_s - 1))
  # Alpha estimado
  if(rho < 0){
    alpha <- 0
  } else{
    alpha <- rho
  }
  # Lambda estimado
  lambda <- (1 - alpha) * y_bar

  # Lista com as informacoes para a saida
  res <- list(call = match.call(),
              method = "kendall",
              order = 1,
              pacf = rho,
              coef = round(alpha, 4),
              lambda = lambda)

  return(res)
}
