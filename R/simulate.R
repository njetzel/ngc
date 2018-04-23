#' Simulate n iid samples from the GGC network with p variables observed over T time points
#' @param n sample size
#' @param edge adjacency matrices
#' @param T number of time points to simulate
#' @param error_sd standard deviation of the error term
#' @param cutt how long to wait for the VAR process to be stationary
#' @return an n x p x T array with the simulated data
#' @export simulate_data
simulate_data <-
  function(
    n,
    edge,
    T,
    error_sd = 1,
    cutt = 30
  ){
    d <- dim(edge)[1]
    p <- dim(edge)[2]
    x = array(0, c(n, p, T+(cutt*d)))

    for (i in 1:n)
    {
      x[i, 1:p, 1:d] <- rnorm(d*p, c(p, d))
      for (j in (d+1): ncol(x[1,,]))
      {
        x[i, ,j] <- rnorm(p, 0, error_sd)
        for (l in 1:d)
        {
          x[i, ,j] <- x[i, ,j] + edge[l,,] %*% x[i, ,j-l]
        }
      }
    }
    return(x[,,-c(1:(cutt*d))])
  }
