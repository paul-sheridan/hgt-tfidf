## multivariate hypergeometric cumulative
pmvhyper <- function(x, n, k, log = FALSE) {
  did_break <- FALSE
  P <- 0
  x1 <- x[1]
  x2 <- x[2]
  n1 <- n[1]
  n2 <- n[2]
  
  suppressWarnings({
  if (x1 == 0 && x2 == 0) {
    x_prime <- c(0, 0, k - 0 - 0)
    P <- P 
  } else if (x1 == 0 && x2 > 1) {
    for (i in 0 : n1) {
      for (j in 0 : (x2 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
      }
      if (did_break == TRUE) { break }
    }
  } else if (x1 > 1 && x2 == 0) {
    for (i in 0 : (x1 - 1)) {
      for (j in 0 : n2) {
        x_prime <- c(i, j, k - i - j)
        #print(paste0("i = ", i, ", j = ", j, ", k = ", k, ", P = ", P, ", log(1-P) = ", log(1-P)))
        P <- P + dmvhyper(x = x_prime, n = n, k = k) 
        if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
      }
      if (did_break == TRUE) { break }
    }
  } else if (x1 > 1 && x2 > 1) {
    for (i in 0 : (x1 - 1)) {
      for (j in 0 : n2) {
        x_prime <- c(i, j, k - i - j)
        P <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
      }
      if (did_break == TRUE) { break }
    }

    for (i in x1 : n1) {
      for (j in 0 : (x2 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
      }
      if (did_break == TRUE) { break }
    }
  } else if (x1 == 0 && x2 == 1) {
    for (i in 0 : n1) {
      x_prime <- c(i, 0, k - i - 0)
      P <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
    }
  } else if (x1 == 1 && x2 == 0) {
    for (j in 0 : n2) {
      x_prime <- c(0, j, k - 0 - j)
      P <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
    }
  } else if (x1 == 1 && x2 == 1) {
    for (i in 0 : n1) {
      x_prime <- c(i, 0, k - i - 0)
      P <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
    }

    for (j in 1 : n2) {
      x_prime <- c(0, j, k - 0 - j)
      P <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
    }
  } else if (x1 == 1 && x2 > 1) {
    for (j in 0 : n2) {
      x_prime <- c(0, j, k - 0 - j)
      P <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
    }

    for (i in 1 : n1) {
      for (j in 0 : (x2 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
      }
      if (did_break == TRUE) { break }
    }
  } else if (x1 > 1 && x2 == 1) {
    for (i in 0 : n1) {
      x_prime <- c(i, 0, k - i - 0)
      P <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
    }

    for (j in 1 : n2) {
      for (i in 0 : (x1 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1-P)) == TRUE) { did_break <- TRUE; break }
      }
      if (did_break == TRUE) { break }
    }
  }

  if (log == FALSE) {
    return(1 - P)
  } else if (log == TRUE && did_break == FALSE) {
    return(-log(1 - P))
  } else if (log == TRUE && did_break == TRUE) {
    return(Inf)
  }
  })
}


## multivariate hypergeometric cumulative
pmvhyper2 <- function(x, n, k, log = FALSE) {
  P <- 0
  x1 <- x[1]
  x2 <- x[2]
  n1 <- n[1]
  n2 <- n[2]
  
  suppressWarnings({
  if (x1 == 0 && x2 == 0) {
    x_prime <- c(0, 0, k - 0 - 0)
    P <- P 
  } else if (x1 == 0 && x2 > 1) {
    for (i in 0 : n1) {
      for (j in 0 : (x2 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1 - P_cand)) == FALSE) { break } else { P <- P_cand }
      }
    }
  } else if (x1 > 1 && x2 == 0) {
    for (i in 0 : (x1 - 1)) {
      for (j in 0 : n2) {
        x_prime <- c(i, j, k - i - j)
        P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
      }
    }
  } else if (x1 > 1 && x2 > 1) {
    for (i in 0 : (x1 - 1)) {
      for (j in 0 : n2) {
        x_prime <- c(i, j, k - i - j)
        P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
      }
    }

    for (i in x1 : n1) {
      for (j in 0 : (x2 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
      }
    }
  } else if (x1 == 0 && x2 == 1) {
    for (i in 0 : n1) {
      x_prime <- c(i, 0, k - i - 0)
      P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
    }
  } else if (x1 == 1 && x2 == 0) {
    for (j in 0 : n2) {
      x_prime <- c(0, j, k - 0 - j)
      P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
    }
  } else if (x1 == 1 && x2 == 1) {
    for (i in 0 : n1) {
      x_prime <- c(i, 0, k - i - 0)
      P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
    }

    for (j in 1 : n2) {
      x_prime <- c(0, j, k - 0 - j)
      P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
    }
  } else if (x1 == 1 && x2 > 1) {
    for (j in 0 : n2) {
      x_prime <- c(0, j, k - 0 - j)
      P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
    }

    for (i in 1 : n1) {
      for (j in 0 : (x2 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
      }
    }
  } else if (x1 > 1 && x2 == 1) {
    for (i in 0 : n1) {
      x_prime <- c(i, 0, k - i - 0)
      P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
      if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
    }

    for (j in 1 : n2) {
      for (i in 0 : (x1 - 1)) {
        x_prime <- c(i, j, k - i - j)
        P_cand <- P + dmvhyper(x = x_prime, n = n, k = k)
        if (is.nan(log(1 - P_cand)) == FALSE) { P <- P_cand } else { break }
      }
    }
  }

  if (log == FALSE) {
    return(1 - P)
  } else if (log == TRUE) {
    return(-log(1 - P))
  }
  })
}





## multivariate hypergeometric cumulative
pmvhyper3 <- function(x, n, k, log = FALSE) {
  P <- 0
  x1 <- x[1]
  x2 <- x[2]
  n1 <- n[1]
  n2 <- n[2]

  if (x1 == 0 && x2 == 0) {
    P <- 1
  } else {
    for (i in x1 : n1) {
      for (j in x2 : n2) {
        x_prime <- c(i, j, k - i - j)
        P <- P + dmvhyper(x = x_prime, n = n, k = k)
      }
    }
  }
  
  if (log == FALSE) {
    return(P)
  } else {
    return(-log(P))
  }
}








