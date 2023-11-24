#------------------------------------------------------------------------------#
# File:        kalman.log.likelihood.R
#
# Description: This function takes as input the coefficient matrices of the
#              given state-space model and the conditional expectation and
#              covariance matrix of the initial state and returns a vector
#              with the log likelihood value from each period,
#              as well as the cumulative sum.
#------------------------------------------------------------------------------#
kalman.log.likelihood <- function(xi.tm1tm1, P.tm1tm1, F, Q, A, H, R, kappa, y, x) {
    t.end <- dim(y)[1]
    n <- dim(y)[2]
    ll.vec <- matrix(NA,t.end,1)
    ll.cum <- 0
    prediction.error.vec <- matrix(NA,t.end,n)
    xi.tt <- xi.tm1tm1
    P.tt  <- P.tm1tm1

    # I have added this line to replace log(2 * atan(1) * 4)
    log_2_pi = 1.837877066409345483560659472811235279722
    # log_2_pi = log(2 * atan(1) * 4)

    # MAKE THE TRANSPOSES OUTSIDE THE LOOP
    t_F = t(F) 
    t_H = t(H)
    t_A = t(A)

    for (t.i in 1:t.end){

        xi.ttm1 <- F %*% xi.tt
        P.ttm1  <- F %*% P.tt %*% t_F + Q
        prediction.error <- (as.vector(y[t.i,]) - as.vector(t_A %*% as.vector(x[t.i,])) - as.vector(t_H %*% xi.ttm1))
        
        HPHR <- t_H %*% P.ttm1 %*% H + (kappa[t.i]^2) * R
        inv_HPHR <- solve(HPHR, prediction.error)

        ll.vec[t.i] <- drop(-(n / 2) * log_2_pi - 0.5 * log(det(HPHR)) -0.5 * prediction.error %*% inv_HPHR)
        ll.cum <- ll.cum + ll.vec[t.i]
        
        kalman.gain <- P.ttm1 %*% H %*% solve(HPHR)
        prediction.error.vec[t.i,] <- prediction.error

        xi.tt <- xi.ttm1 + P.ttm1 %*% H %*% inv_HPHR
        P.tt  <- P.ttm1 -  P.ttm1 %*% H %*% solve(HPHR, t_H %*% P.ttm1)
    }
    return(list("ll.vec"=ll.vec,"ll.cum"=ll.cum,"prediction.error"=prediction.error.vec))
}