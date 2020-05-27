library(expm)

# MC's limit law (if exist), stationary law and mean occupation law
# p is transient probability
MC_limit_law <- function(p) {
    n_row <- nrow(p)
    I <- diag(1, n_row)
    c <- I - p
    c[, n_row] <- matrix(array(1, n_row), nrow = n_row, ncol = 1, byrow = T)
    pi <- solve(c)
    return(pi[n_row, ])
}

# p = transient probability, n = transient step
MC_mean_occupation_time <- function(p, n) {
    x <- p %^% 0
    for (i in c(1:n)) {
        x <- x + p %^% i
    }
    return(x)
}

# p = transient probability, n = target status
MC_hitting_time <- function(p, n) {
    n_row <- nrow(p) - 1
    I <- diag(1, n_row)

    h <- I - p[-n, -n]
    h <- solve(h)
    one <- matrix(array(1, n_row), nrow = n_row, ncol = 1, byrow = T)
    return(h %*% one)
}


# get one-step transition probability of a MP
# Q = transition rate matrix, r = a upper bound of rate
P_hat <- function(Q,r) {
    return(diag(nrow(Q)) + Q / r)
}

# Getting MP's transient probability by using P hat
# p_hat = one-step transition probability, rate = P hat's rate, tm = time, error = the fitting error
transient_probability_by_Phat <- function(p_hat, rate, tm, error) {
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

# Getting MP's transient probability by using Q (rate matrix)
# Q = transition rate matrix, t = time
transition_probability_by_Q <- function(t, Q){
 Pt <- expm(t*Q)
 return(Pt)
 }

# Get MP's mean occupation time
# p_hat =  one-step transition probability, r = p hat's rate, t = time
MP_mean_occupation_time<- function(p_hat,r,t)
{
    Mt<-diag(nrow(p_hat))
    for (n in 1:49) {
        Mt<- Mt+ ppois(n, r*t,lower.tail = FALSE)*p_hat%^%n
    }
    return(r^(-1)*Mt)
}

# j = target status, Q = transition rate martix of MP
MP_hitting_time <- function(j, Q) {
    Qj <- Q[-j, -j]
    inv_Qj <- solve(Qj)
    h_j <- -rowSums(inv_Qj)
    return(h_j)
}

# Q = transition rate martix of MP
MP_limit_law <- function (Q)
{
    n_row = nrow(Q)
    Q[, n_row] = matrix(array(1, n_row), nrow = n_row, ncol = 1, byrow = T)
    Pi <- solve(Q) 
    return(Pi[n_row,])

}

# MOL=mean occupation law or mean occupation time for MP or MC, c=cost which is a col vector
Cost_modle <- function(MOL, c) {
    return(MOL %*% c)
}
