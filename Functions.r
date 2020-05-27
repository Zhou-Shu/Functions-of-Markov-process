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

# estamate MC's transition martrix by giving a sequence of states
# x = the sequence of statuses
MC_est_transition_probs <- function(x) {
    ## let R determine how many states there are,
    ## and replace factor-valued states with integers 1,2,..., nstates
    1
    ## (the following few lines are not necessary,
    ## if your states are already labelled 1, 2, ..., N)
    states <- sort(unique(x))
    nstates <- length(states)
    for (i in 1:nstates)
    {
        x <- replace(x, x == states[i], i)
        x <- as.numeric(x)
    }
    ## Count the number of each type of transition (i -> j),
    ## and represent these counts as a matrix. Note that if
    ## there are length(x) observations, there are only length(x) - 1
    ## transitions to count. (There is no transition
    ## out of the final observation.)
    N <- length(x) - 1
    TransitionCounts <- matrix(0, nrow = nstates, ncol = nstates)
    for (n in 1:N) {
        TransitionCounts[x[n], x[n + 1]] <- TransitionCounts[x[n], x[n + 1]] + 1
    }
    ## Sum the rows of TransitionCounts to give the total number
    ## of transitions out of i. Then divide each row by its
    ## respective total to get the estimated transition probabilities
    P.est <- matrix(0, nrow = nstates, ncol = nstates)
    for (i in 1:nstates) {
        ti <- sum(TransitionCounts[i, ])
        P.est[i, ] <- TransitionCounts[i, ] / ti
    }
    ## Return the estimated transition probability matrix
    return(P.est)
}


# get one-step transition probability of a MP
# Q = generator matrix, r = a upper bound of rate
P_hat <- function(Q,r) {
    return(diag(nrow(Q)) + Q / r)
}

# Getting MP's transient probability by using P hat
# p_hat = one-step transition probability, rate = P hat's rate, tm = time, error = the fitting error
transition_probability_by_Phat <- function(p_hat, rate, tm, error) {
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

# Getting MP's transient probability by using P hat
# p_hat = one-step transition probability, rate = P hat's rate, tm = time, M = the number of iterations
transition_probability_by_Phat_no_error <- function(p_hat, rate, tm, M) {
    nstates <- nrow(p_hat)
    Pt <- diag(nstates)
    p_hat_n <- p_hat
    mean_jumps <- rate * tm
    for (n in 1:M)
    {
        Pt <- Pt + ((mean_jumps)^n / factorial(n)) * p_hat_n
        p_hat_n <- p_hat_n %*% p_hat
    }
    Pt <- exp(-mean_jumps) * Pt
    return(Pt)
}

# Getting MP's transient probability by using Q (rate matrix)
# Q = generator matrix, t = time
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

# j = target status, Q = generator martix of MP
MP_hitting_time <- function(j, Q) {
    Qj <- Q[-j, -j]
    inv_Qj <- solve(Qj)
    h_j <- -rowSums(inv_Qj)
    return(h_j)
}

# Q = generator martix of MP
MP_limit_law <- function (Q)
{
    n_row = nrow(Q)
    Q[, n_row] = matrix(array(1, n_row), nrow = n_row, ncol = 1, byrow = T)
    Pi <- solve(Q) 
    return(Pi[n_row,])

}

# Simulate data for the Markov process in the form (xi,ti-t_{i-1} ,ti) for i = 0, 1, . . . , 200, where xi is the MP visits in the ith jump and ti is the time that MP stays in state xi before the next jump to another state.
# Q = generator matrix, init_state = initial state, N = # of step
MP_Simulate_data <- function(Q, init_state, N) {
    ## Define the jump rates
    r <- c()
    for (i in c(1:nrow(Q))) {
        r <- cbind(r, -Q[i, i])
    }
    P <- P_hat(Q, max(r))
    ## Set the initial values
    state <- init_state
    t <- 0
    ## Visits to any state have duration Exp(r_i), and this hold time
    ## is followed by departure to the next state, which has
    ## conditional distribution P[i,] given the present state i.
    ## The MP sample path is stored as a 3-row matrix, with the first
    ## row giving the sequence of visited states, the second row giving
    ## the hold times for each visit, and the third row giving
    ## cumulative time.
    MP <- matrix(c(state, 0, t), ncol = 1)
    for (i in 1:N) {
        hold.time <- rexp(1, rate = r[state])
        t <- t + hold.time
        MP <- cbind(MP, c(state, hold.time, t))
        state <- sample(c(1:nrow(Q)), size = 1, prob = P[state, ])
    }
    return(MP)
}

# MOL=mean occupation law or mean occupation time for MP or MC, c=cost which is a col vector
Cost_modle <- function(MOL, c) {
    return(MOL %*% c)
}
