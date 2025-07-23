#' @title Poisson INAR(1) Forecast
#'
#' @rdname forecast_inar1
#'
#' @param series vector of (non-negative) integer time series
#' @param est named list of (non-negative) estimated values of alpha and lambda.
#'
#' @importFrom stats dbinom dpois
#' @references Freeland, R. K. & McCabe, Brendan (2004).
#' Forecasting discrete valued low count times series.
#'
#' @examples
#' # Estimated values with Yule-Walker method
#' sim <- sim_pois(100, 0.3, 1)
#' est <- yw(sim)
#' forecast_inar(sim, est)
#' @export
forecast_inar1 <- function(series, est){
  # Probabilidade Acumulada
  prob_forec_inar <- function(x, x_t, est_alpha, est_lambda){
    # x -> valor testado pra previsao; x_t -> ultimo valor da serie temporal
    cm <- 0
    for(s in 0:min(x, x_t)){
      prob <- dbinom(s, x_t, est_alpha) * dpois((x - s), est_lambda)
      cm <- cm + prob
    }

    return(cm)
  }

  prob_cm <- 0 # probabilidade acumulada
  x <- 0 # valor inicial para a predicao
  while(prob_cm <= 0.5){ # calculando a probabilidade acumulada ate a mediana
    prob_cm <- prob_cm + prob_forec_inar(x, series[length(series)],
                                         est$coef, est$lambda)
    x <- x + 1 # valor inteiro
  }

  return(x - 1) # a probabilidade acumulada provavelmente passa de 0.5
}
