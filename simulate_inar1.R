#' @title Simulation from a Poisson INAR(1)
#'
#' @rdname simulate_inar1
#'
#' @param n number of observations. If length(n) > 1, the length is taken to
#' be the number required.
#' @param alpha vector of probabilities
#' @param lambda vector of (non-negative) means.
#' @param n.start length of ‘burn-in’ period. If NA, the default, a reasonable
#' value is computed.
#'
#' @return `sim_pois` returns an time series object of class "`ts`".
#'
#' @importFrom stats ts rbinom rpois
#' @references Alzaid, A. A., & Al-Osh, M. A. (1987).
#' First-order integer-valued autoregressive (INAR (1)) process.
#'
#' @examples
#' sim_pois(100, 0.3, 1)
#'
#' plot(sim_pois(300, 0.4, 2))
#' @export
sim_pois <- function(n, alpha, lambda, n.start = 100) { # binomial thinning

  if(n%%1 != 0) stop("n is not (non-negative) integer")
  if(!is.numeric(alpha) ) stop("alpha is not numeric")
  if(!is.numeric(lambda)) stop("lambda is not numeric")

  # Tamanho da serie + burnin
  len <- n + n.start

  # Erros com distribuicao Poisson e lambda definido
  error <- rpois(len, lambda)

  # Vetor do tamanho da serie simulada + burnin
  sim <- rep(0, times = len) # primeiro valor = primeiro erro
  sim[1] <- error[1] # por causa do burnin nao tem importancia

  for(i in 2:len){ # Equacao 2.1 ; M. A. AL-OSH AND A. A. ALZAI (1987)
    sim[i] <- rbinom(1, sim[i-1], alpha) + error[i]
  }

  sim <- ts(sim[(n.start + 1):len])

  return(sim)
}
