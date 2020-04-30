library(expm)

limit_law <- function(p) {
    n_row <- nrow(p)
    I <- diag(1, n_row)
    c <- I - p
    c[, n_row] <- matrix(array(1, n_row), nrow = n_row, ncol = 1, byrow = T)
    pi <- solve(c)
    return(pi[n_row, ])
}


mean_occupation_time <- function(p, n) {
    x <- p %^% 0
    for (i in c(1:n)) {
        x <- x + p %^% i
    }
    return(x)
}

hitting_time <- function(p, n) {
    n_row <- nrow(p) - 1
    I <- diag(1, n_row)

    h <- I - p[-n, -n]
    h <- solve(h)
    one <- matrix(array(1, n_row), nrow = n_row, ncol = 1, byrow = T)
    return(h %*% one)
}

transient_probability <- function(p_hat, rate, tm, error) {
    nstates <- nrow(p_hat)
    Pt <- diag(nstates)
    p_hat_n <- p_hat
    mean_jumps <- rate * tm
    M <- qpois(error, mean_jumps, lower.tail = FALSE)
    for (n in 1:M)
    {
        Pt <- Pt + ((mean_jumps)^n / factorial(n)) * p_hat_n
        p_hat_n <- p_hat_n %*% p_hat
    }
    Pt <- exp(-mean_jumps) * Pt
    return(Pt)
}


MP_mean_occupation_time<- function(p_hat,r,t)
{
    Mt<-diag(nrow(p_hat))
    for (n in 1:49) {
        Mt<- Mt+ ppois(n, r*t,lower.tail = FALSE)*p_hat%^%n
    }
    return(r^(-1)*Mt)
}
