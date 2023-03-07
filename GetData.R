GenData <- function(nprim = 200,
                    naux = 500, 
                    setting = 1) {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    X1_a = runif(naux, 0, 1)
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    b1_p = -0.5
    b2_p = 0.5
    XB_p = cbind(X1_p, X2_p) %*% c(b1_p, b2_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 0.5
    } else if(setting %in% c(2,4)) {
        b1_a = -0.5
        b2_a = 0.2
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)*3/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)

    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                          X2 = X2_a,
                          time = survTime_a,
                          status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}

GenData_LargerEff <- function(nprim = 200,
                    naux = 500, 
                    setting = 1) {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    X1_a = runif(naux, 0, 1)
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    b1_p = -0.5
    b2_p = 2
    XB_p = cbind(X1_p, X2_p) %*% c(b1_p, b2_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 2
    } else if(setting %in% c(2,4)) {
        b1_a = -0.5
        b2_a = 0.2
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)*5/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)
    
    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                         X2 = X2_a,
                         time = survTime_a,
                         status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}

getCumEval <- function(bhest) {
    breakPoint <- seq(0.15, 2.5, length.out = 50)
    Heval <- rep(NA, length(breakPoint))
    for(i in 1:length(breakPoint)) {
        diff <- breakPoint[i] - bhest$time
        Heval[i] <- bhest$hazard[which.min(abs(diff))]
    }
    return(Heval)
}


GenData_DiffDist <- function(nprim = 200,
                    naux = 500, 
                    setting = 1) {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    X1_a = runif(naux, -0.5, 0.5)
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    b1_p = -0.5
    b2_p = 0.5
    XB_p = cbind(X1_p, X2_p) %*% c(b1_p, b2_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 0.5
    } else if(setting %in% c(2,4)) {
        b1_a = -0.5
        b2_a = 0.2
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)*3/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)
    
    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                         X2 = X2_a,
                         time = survTime_a,
                         status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}




GenData_DiffDist2 <- function(nprim = 200,
                             naux = 500, 
                             setting = 1,
                             Dist = "Norm") {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    if (Dist == "Norm") {
        X1_a = rnorm(naux, 0, 1)
    } else if (Dist == "Beta") {
        X1_a = rbeta(naux, 1, 2)
    }
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    b1_p = -0.5
    b2_p = 0.5
    XB_p = cbind(X1_p, X2_p) %*% c(b1_p, b2_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 0.5
    } else if(setting %in% c(2,4)) {
        b1_a = -0.5
        b2_a = 0.2
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a) %*% c(b1_a, b2_a)
        scale_a = exp(-XB_a/2)*3/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)
    
    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                         X2 = X2_a,
                         time = survTime_a,
                         status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}


GenData_moreCov <- function(nprim = 200,
                    naux = 500, 
                    setting = 1,
                    Dist = "Norm") {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    X1_a = runif(naux, 0, 1)
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    if (Dist == "Norm") {
        X3_p = rnorm(nprim, 0, 1)
        X3_a = rnorm(naux, 0.5, 1)
    } else if (Dist == "Beta") {
        X3_p = runif(nprim, 0, 1)
        X3_a = rbeta(naux, 1, 2)
    }

    b1_p = -0.5
    b2_p = 0.5
    b3_p = 0.1
    
    XB_p = cbind(X1_p, X2_p, X3_p) %*% c(b1_p, b2_p, b3_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 0.5
        b3_a = 0.1
    } else if(setting %in% c(2,4)) {
        b1_a = -0.5
        b2_a = 0.2
        b3_a = 0.1
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a, X3_a) %*% c(b1_a, b2_a, b3_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a, X3_a) %*% c(b1_a, b2_a, b3_a)
        scale_a = exp(-XB_a/2)*3/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)
    
    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          X3 = X3_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                         X2 = X2_a,
                         X3 = X3_a,
                         time = survTime_a,
                         status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}




GenData_FiveCov <- function(nprim = 250,
                            naux = 6400, 
                            setting = 1,
                            Dist = "Norm") {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    X1_a = runif(naux, 0, 1)
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    if (Dist == "Norm") {
        X3_p = rnorm(nprim, 0, 1)
        X3_a = rnorm(naux, 0.5, 1)
    } else if (Dist == "Beta") {
        X3_p = runif(nprim, 0, 1)
        X3_a = rbeta(naux, 1, 2)
    }
    
    X4_p = runif(nprim, 0, 1)
    X4_a = runif(naux, 0, 1)
    
    X5_p = runif(nprim, 0, 1)
    X5_a = runif(naux, 0, 1)
    
    b1_p = -0.5
    b2_p = 0.5
    b3_p = 0.2
    b4_p = 0.1
    b5_p = 0.1
    
    XB_p = cbind(X1_p, X2_p, X3_p, X4_p, X5_p) %*% c(b1_p, b2_p, b3_p, b4_p, b5_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 0.5
        b3_a = 0.2
        b4_a = 0.1
        b5_a = 0.1
    } else if(setting %in% c(2,4)) {
        b1_a = -0.5
        b2_a = 0.2
        b3_a = 0.2
        b4_a = 0.1
        b5_a = 0.1
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a, X3_a, X4_a, X5_a) %*% c(b1_a, b2_a, b3_a, b4_a, b5_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a, X3_a, X4_a, X5_a) %*% c(b1_a, b2_a, b3_a, b4_a, b5_a)
        scale_a = exp(-XB_a/2)*3/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)
    
    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          X3 = X3_p,
                          X4 = X4_p,
                          X5 = X5_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                         X2 = X2_a,
                         X3 = X3_a,
                         X4 = X4_a,
                         X5 = X5_a,
                         time = survTime_a,
                         status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}


GenData_FiveCov_NonSparse <- function(nprim = 250,
                            naux = 6400, 
                            setting = 1,
                            Dist = "Norm") {
    
    ### two covarites first
    X1_p = runif(nprim, 0, 1)
    X2_p  = rbinom(nprim, size = 1, prob = 0.5)
    
    X1_a = runif(naux, 0, 1)
    X2_a = rbinom(naux, size = 1, prob = 0.5)
    
    if (Dist == "Norm") {
        X3_p = rnorm(nprim, 0, 1)
        X3_a = rnorm(naux, 0.5, 1)
    } else if (Dist == "Beta") {
        X3_p = runif(nprim, 0, 1)
        X3_a = rbeta(naux, 1, 2)
    }
    
    X4_p = runif(nprim, 0, 1)
    X4_a = runif(naux, 0, 1)
    
    X5_p = runif(nprim, 0, 1)
    X5_a = runif(naux, 0, 1)
    
    b1_p = -0.5
    b2_p = 0.5
    b3_p = 0.2
    b4_p = 0.1
    b5_p = 0.1
    
    XB_p = cbind(X1_p, X2_p, X3_p, X4_p, X5_p) %*% c(b1_p, b2_p, b3_p, b4_p, b5_p)
    scale_p = exp(-XB_p/2)
    shape_p = 2
    
    if(setting %in% c(1,3)) {
        b1_a = -0.5
        b2_a = 0.5
        b3_a = 0.2
        b4_a = 0.1
        b5_a = 0.1
    } else if(setting %in% c(2,4)) {
        b1_a = -0.7
        b2_a = 0.2
        b3_a = 0.2
        b4_a = -0.1
        b5_a = 0.1
    }
    
    if(setting %in% c(1,2)) {
        XB_a = cbind(X1_a, X2_a, X3_a, X4_a, X5_a) %*% c(b1_a, b2_a, b3_a, b4_a, b5_a)
        scale_a = exp(-XB_a/2)
        shape_a = 2
    } else if(setting %in% c(3,4)) {
        XB_a = cbind(X1_a, X2_a, X3_a, X4_a, X5_a) %*% c(b1_a, b2_a, b3_a, b4_a, b5_a)
        scale_a = exp(-XB_a/2)*3/2
        shape_a = 2
    }
    
    T_p = rweibull(nprim, shape = shape_p, scale = scale_p)
    C_p = runif(nprim, 0, 4.55)
    survTime_p = ifelse(T_p < C_p, T_p, C_p)
    event_p = ifelse(T_p < C_p, 2, 1)
    
    T_a = rweibull(naux, shape = shape_a, scale = scale_a)
    C_a = runif(naux, 0, 4.55)
    survTime_a = ifelse(T_a < C_a, T_a, C_a)
    event_a = ifelse(T_a < C_a, 2, 1)
    
    primData = data.frame(X1 = X1_p,
                          X2 = X2_p,
                          X3 = X3_p,
                          X4 = X4_p,
                          X5 = X5_p,
                          time = survTime_p,
                          status = event_p)
    auxData = data.frame(X1 = X1_a,
                         X2 = X2_a,
                         X3 = X3_a,
                         X4 = X4_a,
                         X5 = X5_a,
                         time = survTime_a,
                         status = event_a)
    
    return(list(primData = primData,
                auxData = auxData))
}





