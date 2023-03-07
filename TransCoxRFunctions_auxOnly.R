dQtocumQ <- function(dQ, status = NULL) {
    if(length(status)>0) {
        cumQ = rep(NA, length(status))
        ncount = 1
        for(i in 1:length(status)) {
            if(status[i] == 2) {
                cumQ[i] = sum(dQ[1:ncount])
                ncount = ncount + 1
                if(ncount > length(dQ)) {
                    ncount = length(dQ)
                }
            } else {
                cumQ[i] = sum(dQ[1:ncount])
            }
        }
    } else {
        cumQ = rep(NA, length(dQ))
        for(i in 1:length(dQ)) {
            cumQ[i] = sum(dQ[1:i])
        }
    }
    return(cumQ)
}
runTransCox_one <- function(Pout, l1 = 1, l2 = 1, learning_rate = 0.004, nsteps = 200,
                            cov = c('X1', 'X2')){
    CovData = Pout$primData[, cov]
    status = Pout$primData[, "status"]
    cumH = Pout$primData$fullCumQ
    hazards = Pout$dQ$dQ
    
    test <- TransCox(CovData = as.matrix(CovData), 
                     cumH = cumH, 
                     hazards = hazards, 
                     status = status, 
                     estR = Pout$estR, 
                     Xinn = Pout$Xinn, 
                     lambda1 = l1, lambda2 = l2,
                     learning_rate = learning_rate,
                     nsteps = nsteps)
    names(test) <- c("eta", "xi")
    
    return(list(eta = test$eta,
                xi = test$xi,
                new_beta = Pout$estR + test$eta,
                new_IntH = Pout$dQ$dQ + test$xi,
                time = Pout$primData[status == 2, "time"]))
}

GetLogLike <- function(status, CovData, hazards, 
                       newBeta,
                       newHaz) {
    XB = as.matrix(CovData) %*% newBeta
    newCumH <- dQtocumQ(newHaz, status)
    LogL <- sum(XB[status == 2] + log(newHaz)) - sum(newCumH * exp(XB))
    return(LogL)
}

GetBIC <- function(status, CovData, hazards, 
                   newBeta,
                   newHaz,
                   eta,
                   xi,cutoff = 1e-5) {
    Logl <- GetLogLike(status = status, 
                       CovData = CovData,
                       hazards = hazards,
                       newBeta = newBeta,
                       newHaz = newHaz)
    K = sum(c(abs(eta), abs(xi))>cutoff)
    n = length(status)
    BIC <- K*log(n) - 2*Logl
    return(BIC)
}

CheckLR_BIC <- function(primData, auxData, cov = c("X1", "X2"),
                        statusvar = "status",
                        lr_vec = seq(0.001, 0.004, by = 0.0005),
                        l1 = 0.1, l2 = 0.1,
                        nsteps = 200) {
    Cout <- GetAuxSurv(auxData, cov = cov)
    Pout <- GetPrimaryParam(primData, q = Cout$q, estR = Cout$estR)
    
    CovData = Pout$primData[, cov]
    status = Pout$primData[, statusvar]
    cumH = Pout$primData$fullCumQ
    hazards = Pout$dQ$dQ
    
    BICvec <- rep(NA, length(lr_vec))
    for(i in 1:length(lr_vec)) {
        test <- TransCox(CovData = as.matrix(CovData), 
                         cumH = cumH, 
                         hazards = hazards, 
                         status = status, 
                         estR = Pout$estR, 
                         Xinn = Pout$Xinn, 
                         lambda1 = l1, lambda2 = l2,
                         learning_rate = lr_vec[i],
                         nsteps = nsteps)
        names(test) <- c("eta", "xi")
        newBeta = Pout$estR + test$eta
        newHaz = Pout$dQ$dQ + test$xi
        
        BICvalue <- GetBIC(status = status, 
                           CovData = CovData,
                           hazards = hazards,
                           newBeta = newBeta,
                           newHaz = newHaz,
                           eta = test$eta,
                           xi = test$xi,
                           cutoff = 1e-5)
        BICvec[i] <- BICvalue
    }
    names(BICvec) = lr_vec
    
    idx0 <- which.min(BICvec)
    
    b_lr <- lr_vec[idx0]
    
    return(list(best_lr = b_lr,
                BICvec = BICvec))
}
SelParam_By_BIC <- function(primData, auxData, cov = c("X1", "X2"),
                            statusvar = "status",
                            lambda1_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                            lambda2_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                            learning_rate = 0.004,
                            nsteps = 100) {
    Cout <- GetAuxSurv(auxData, cov = cov)
    Pout <- GetPrimaryParam(primData, q = Cout$q, estR = Cout$estR)
    
    CovData = Pout$primData[, cov]
    status = Pout$primData[, statusvar]
    cumH = Pout$primData$fullCumQ
    hazards = Pout$dQ$dQ
    
    BICmat <- matrix(NA, length(lambda1_vec), length(lambda2_vec))
    pb = txtProgressBar(min = 0, max = (length(lambda1_vec))^2, style = 2, initial = 0) 
    for(i in 1:length(lambda1_vec)) {
        for(j in 1:length(lambda2_vec)) {
            
            thisi = (i-1) * length(lambda1_vec) + j
            setTxtProgressBar(pb, thisi)
            
            lambda1 = lambda1_vec[i]
            lambda2 = lambda2_vec[j]
            
            test <- TransCox(CovData = as.matrix(CovData), 
                             cumH = cumH, 
                             hazards = hazards, 
                             status = status, 
                             estR = Pout$estR, 
                             Xinn = Pout$Xinn, 
                             lambda1 = lambda1, lambda2 = lambda2,
                             learning_rate = learning_rate,
                             nsteps = nsteps)
            names(test) <- c("eta", "xi")
            newBeta = Pout$estR + test$eta
            newHaz = Pout$dQ$dQ + test$xi
            
            BICvalue <- GetBIC(status = status, 
                               CovData = CovData,
                               hazards = hazards,
                               newBeta = newBeta,
                               newHaz = newHaz,
                               eta = test$eta,
                               xi = test$xi,
                               cutoff = 1e-5)
            
            BICmat[i,j] <- BICvalue
        }
    }
    close(pb)
    rownames(BICmat) = lambda1_vec
    colnames(BICmat) <- lambda2_vec
    
    idx0 <- which(BICmat == min(BICmat, na.rm = TRUE), arr.ind = TRUE)
    
    b_lambda1 <- lambda1_vec[idx0[1]]
    b_lambda2 <- lambda2_vec[idx0[2]]
    
    return(list(best_la1 = b_lambda1,
                best_la2 = b_lambda2,
                BICmat = BICmat))
}



SelLR_By_BIC <- function(primData, auxData, cov = c("X1", "X2"),
                            statusvar = "status",
                            lambda1 = 0.1,
                            lambda2 = 0.1,
                            learning_rate_vec = c(0.0001, 0.0005, 0.001, 0.002,0.003, 0.004, 0.005),
                            nsteps_vec = c(100,200,300,400)) {
    Cout <- GetAuxSurv(auxData, cov = cov)
    Pout <- GetPrimaryParam(primData, q = Cout$q, estR = Cout$estR)
    
    CovData = Pout$primData[, cov]
    status = Pout$primData[, statusvar]
    cumH = Pout$primData$fullCumQ
    hazards = Pout$dQ$dQ
    
    BICmat <- matrix(NA, length(learning_rate_vec), length(nsteps_vec))
    pb = txtProgressBar(min = 0, max = (length(learning_rate_vec))^2, style = 2, initial = 0) 
    for(i in 1:length(learning_rate_vec)) {
        for(j in 1:length(nsteps_vec)) {
            
            thisi = (i-1) * length(learning_rate_vec) + j
            setTxtProgressBar(pb, thisi)
            
            learning_rate = learning_rate_vec[i]
            nsteps = nsteps_vec[j]
            
            test <- TransCox(CovData = as.matrix(CovData), 
                             cumH = cumH, 
                             hazards = hazards, 
                             status = status, 
                             estR = Pout$estR, 
                             Xinn = Pout$Xinn, 
                             lambda1 = lambda1, 
                             lambda2 = lambda2,
                             learning_rate = learning_rate,
                             nsteps = nsteps)
            names(test) <- c("eta", "xi")
            newBeta = Pout$estR + test$eta
            newHaz = Pout$dQ$dQ + test$xi
            
            BICvalue <- GetBIC(status = status, 
                               CovData = CovData,
                               hazards = hazards,
                               newBeta = newBeta,
                               newHaz = newHaz,
                               eta = test$eta,
                               xi = test$xi,
                               cutoff = 1e-5)
            
            BICmat[i,j] <- BICvalue
        }
    }
    close(pb)
    rownames(BICmat) = learning_rate_vec
    colnames(BICmat) <- nsteps_vec
    
    idx0 <- which(BICmat == min(BICmat, na.rm = TRUE), arr.ind = TRUE)
    
    b_learning_rate <- learning_rate_vec[idx0[1]]
    b_nsteps <- nsteps_vec[idx0[2]]
    
    return(list(best_lr = b_learning_rate,
                best_nsteps = b_nsteps,
                BICmat = BICmat))
}


deltaQ <- function(primData, q) {
    #### get cumulative hazards for primary data
    #### from the combined data
    obsData <- primData[primData$status == 2, ]
    obsData <- obsData[order(obsData$time), ]
    newQ <- data.frame(dQ  = rep(NA, nrow(obsData)),
                       cumQ = rep(NA, nrow(obsData)),
                       t = rep(NA, nrow(obsData)))
    for(i in 1:nrow(obsData)) {
        if (i == 1) {
            newQ$t[i] <- obsData$time[i]
            tmp1 <- q$breakPoints<=obsData$time[i]
            if(all(!tmp1)) {
                idx0 <- 1
            } else {
                idx0 <- max(which(tmp1))
            }
            newQ$dQ[i] = newQ$cumQ[i] <- q$cumHazards[idx0]
        } else {
            newQ$t[i] <- obsData$time[i]
            tmp1 <- q$breakPoints<=obsData$time[i]
            if(all(!tmp1)) {
                idx0 <- 1
            } else {
                idx0 <- max(which(tmp1))
            }
            newQ$cumQ[i] <- q$cumHazards[idx0]
            newQ$dQ[i] <- q$cumHazards[idx0] - newQ$cumQ[(i-1)]
        }
    }
    newQ$dQ[newQ$dQ == 0] <- 0.0001
    newQ$cumQ_upd <- newQ$cumQ
    for(i in 1:nrow(newQ)) {
        newQ$cumQ_upd[i] <- sum(newQ$dQ[1:i])
    }
    return(newQ)
}

GetAuxSurv <- function(auxData, weights = NULL, cov = c("X1", "X2")) {
    functext = paste0("res.cox <- coxph(Surv(time, status == 2) ~ ", paste(cov, collapse = "+"), ", data = auxData, weights = weights)")
    ### Give auxiliary data, estimate r and q
    eval(parse(text = functext))
    bhest <- basehaz(res.cox, centered=FALSE) ## get baseline cumulative hazards
    estR <- res.cox$coefficients
    q <- data.frame(cumHazards = bhest$hazard,
                    breakPoints = bhest$time)
    return(list(estR = estR,
                q = q))
}

GetPrimaryParam <- function(primData, q, estR) {
    ### get primary data-specific breakpoints and hazards
    primData <- primData[order(primData$time), ]
    dQ <- deltaQ(primData, q)
    primData <- primData[order(primData$time), ]
    fullcum <- rep(0, nrow(primData))
    idx0 <- match(dQ$t, primData$time)
    
    ### get cumQ
    fullcum[idx0] <- dQ$cumQ_upd
    for(i in 1:length(fullcum)) {
        if(fullcum[i] == 0 & i>1) {
            fullcum[i] = fullcum[i-1]
        }
    }
    primData$fullCumQ = fullcum

    #### get Ximat
    Ximat <- matrix(0, nrow(primData), nrow(primData))
    for(i in 1:nrow(primData)) {
        tmpidx = rep(0, nrow(primData))
        for(j in 1:i) {
            if(primData$status[j] == 2)
                tmpidx[j] = 1
        }
        Ximat[i,] <- tmpidx
    }
    Xinn <- Ximat[,which(primData$status == 2)]
    return(list(primData = primData,
                Xinn = Xinn,
                dQ = dQ,
                estR = estR))
}

GetSurvEval <- function(pData, myCum, myB, evalT = c(1,2)) {
    # medianCov = c(median(pData$X1), as.numeric(names(which.max(table(pData$X2)))))
    medianCov = c(0.5, 1)
    newCum = data.frame(H = myCum$hazard * c(exp(medianCov %*% myB)),
                        T = myCum$time)
    newCum$SurvF = exp(-newCum$H)
    newCum <- newCum[order(newCum$T), ]
    newCum <- rbind(c(0,0,1), newCum)
    if(max(newCum$T)<2) {
        newCum <- rbind(newCum, c(max(newCum$H),2,0))
    }
    
    ## get RMST
    func = approxfun(newCum$T, newCum$SurvF)
    RMST = rep(NA, length(evalT))
    for(i in 1:length(evalT)) {
        RMST[i] = integrate(func, 0, evalT[i], subdivisions=2000, stop.on.error = FALSE)$value
    }
    return(list(newCum = newCum,
                RMST = RMST))
}

GetSurvEval_3Cov <- function(pData, myCum, myB, evalT = c(1,2)) {
    # medianCov = c(median(pData$X1), as.numeric(names(which.max(table(pData$X2)))))
    medianCov = c(0.5, 1, 0.1)
    newCum = data.frame(H = myCum$hazard * c(exp(medianCov %*% myB)),
                        T = myCum$time)
    newCum$SurvF = exp(-newCum$H)
    newCum <- newCum[order(newCum$T), ]
    newCum <- rbind(c(0,0,1), newCum)
    if(max(newCum$T)<2) {
        newCum <- rbind(newCum, c(max(newCum$H),2,0))
    }
    
    ## get RMST
    func = approxfun(newCum$T, newCum$SurvF)
    RMST = rep(NA, length(evalT))
    for(i in 1:length(evalT)) {
        RMST[i] = integrate(func, 0, evalT[i], subdivisions=2000, stop.on.error = FALSE)$value
    }
    return(list(newCum = newCum,
                RMST = RMST))
}


GetSurvEval_5Cov <- function(pData, myCum, myB, evalT = c(1,2)) {
    # medianCov = c(median(pData$X1), as.numeric(names(which.max(table(pData$X2)))))
    medianCov = c(0.5, 1, 0.1, 0.5, 0.5)
    newCum = data.frame(H = myCum$hazard * c(exp(medianCov %*% myB)),
                        T = myCum$time)
    newCum$SurvF = exp(-newCum$H)
    newCum <- newCum[order(newCum$T), ]
    newCum <- rbind(c(0,0,1), newCum)
    if(max(newCum$T)<2) {
        newCum <- rbind(newCum, c(max(newCum$H),2,0))
    }
    
    ## get RMST
    func = approxfun(newCum$T, newCum$SurvF)
    RMST = rep(NA, length(evalT))
    for(i in 1:length(evalT)) {
        RMST[i] = integrate(func, 0, evalT[i], subdivisions=2000, stop.on.error = FALSE)$value
    }
    return(list(newCum = newCum,
                RMST = RMST))
}


runPerm_one <- function(Pout, weight, l1 = 1, l2 = 1, learning_rate = 0.004, nsteps = 200){
    CovData = Pout$primData[, c('X1', 'X2')]
    status = Pout$primData[, "status"]
    cumH = Pout$primData$fullCumQ
    hazards = Pout$dQ$dQ
    
    test <- TransCox_perm(CovData = as.matrix(CovData),
                          weight = weight,
                     cumH = cumH, 
                     hazards = hazards, 
                     status = status, 
                     estR = Pout$estR, 
                     Xinn = Pout$Xinn, 
                     lambda1 = l1, lambda2 = l2,
                     learning_rate = learning_rate,
                     nsteps = nsteps)
    names(test) <- c("eta", "xi")
    
    return(list(eta = test$eta,
                xi = test$xi,
                new_beta = Pout$estR + test$eta,
                new_IntH = Pout$dQ$dQ + test$xi,
                time = Pout$primData[status == 2, "time"]))
}

runPerm_full <- function(aData,
                         pData,
                         twopoint,
                         npermutation = 100) {
    rec_beta = rec_twop_cumH = rec_twop_surv <- matrix(NA, npermutation, 2)
    rec_cumH <- matrix(NA, npermutation, 50)
    rec_survTrans = matrix(NA, npermutation, 50)
    
    pb = txtProgressBar(min = 0, max = npermutation, style = 2, initial = 0) 
    for(i in 1:npermutation) {
        setTxtProgressBar(pb, i)
        
        aweights <- rexp(nrow(aData), rate = 1)
        Cout <- GetAuxSurv(aData, weights = aweights)
        Pout <- GetPrimaryParam(pData, q = Cout$q, estR = Cout$estR)
        pweights <- rexp(dim(pData)[1], rate = 1)
        Tres_btsp <- runPerm_one(Pout, weight = pweights, 
                                 l1 = BestLam1, l2 = BestLam2, 
                                 learning_rate = learning_rate, nsteps = nsteps)

        Htime <- Tres_btsp$time
        DupN <- sum(duplicated(Htime)) 
        Htime0 <- Tres_btsp$time
        while (DupN >= 1){
            Htime[duplicated(Htime0)] <- Htime0[duplicated(Htime0)] + 0.0001
            DupN <- sum(duplicated(Htime)) 
            Htime0 = Htime
        }
        myHvec <- data.frame(hazard = dQtocumQ(Tres_btsp$new_IntH), time = Htime)
        t1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = twopoint)$y
        ttt1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = seq(0.15, 2, length.out = 50))$y
        
        point50 = seq(0.15, 2, length.out = 50)^2
        bhest_true = data.frame(time = point50, 
                                hazards = point50^2)
        TransSurvRes <- GetSurvEval(pData = pData, myCum = myHvec, myB = Tres$new_beta, evalT = twopoint)
        SurvRes_true <- GetSurvEval(pData = pData, myCum = bhest_true, myB = trueBeta, evalT = twopoint)
        surv1 <- approx(x = TransSurvRes$newCum$T, y = TransSurvRes$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        survTrue <- approx(x = SurvRes_true$newCum$T, y = SurvRes_true$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        
        rec_beta[i,] <- Tres_btsp$new_beta
        rec_twop_cumH[i,] <- t1
        rec_cumH[i,] <- ttt1
        rec_twop_surv[i,] <- TransSurvRes$RMST
        rec_survTrans[i,] <- surv1
    }
    close(pb)
    return(list(rec_beta = rec_beta,
                beta_SE_btsp = apply(rec_beta, 2, sd),
                rec_twop_cumH = rec_twop_cumH,
                cumH_twop_SE_btsp = apply(rec_twop_cumH, 2, sd),
                rec_cumH = rec_cumH,
                cumH_SE_btsp = apply(rec_cumH, 2, sd),
                rec_twop_surv = rec_twop_surv,
                surv_twop_SE_btsp = apply(rec_twop_surv, 2, sd),
                rec_survTrans = rec_survTrans,
                surv_SE_btsp = apply(rec_survTrans, 2, sd)))
}


runBtsp_transCox <- function(aData,
                             pData,
                             twopoint,
                             Tres,
                             BestLam1 = 0.1,
                             BestLam2 = 0.1,
                    nbootstrap = 100,
                    trueBeta = c(-0.5, 0.5),
                    cov = c("X1", "X2")) {
    
    rec_beta = rec_twop_cumH = rec_twop_surv <- matrix(NA, nbootstrap, 2)
    rec_cumH <- matrix(NA, nbootstrap, 50)
    rec_survTrans = matrix(NA, nbootstrap, 50)
    
    pb = txtProgressBar(min = 0, max = nbootstrap, style = 2, initial = 0) 
    for(i in 1:nbootstrap) {
        setTxtProgressBar(pb, i)
        
        aData_btsp <- aData[sample(1:nrow(aData), nrow(aData), replace = TRUE), ]
        pData_btsp <- pData[sample(1:nrow(pData), nrow(pData), replace = TRUE), ]
        Cout_btsp <- GetAuxSurv(aData_btsp, cov = cov)
        Pout_btsp <- GetPrimaryParam(pData_btsp, q = Cout_btsp$q, estR = Cout_btsp$estR)
        Tres_btsp <- runTransCox_one(Pout_btsp, 
                                     l1 = BestLam1, 
                                     l2 = BestLam2, 
                                     learning_rate = learning_rate, 
                                     nsteps = nsteps,
                                     cov = cov)
        Htime <- Tres_btsp$time
        DupN <- sum(duplicated(Htime)) 
        myHvec <- data.frame(hazard = dQtocumQ(Tres_btsp$new_IntH), time = Htime)
        t1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = twopoint)$y
        ttt1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = seq(0.15, 2, length.out = 50))$y
        
        point50 = seq(0.15, 2, length.out = 50)^2
        bhest_true = data.frame(time = point50, 
                                hazards = point50^2)
        # TransSurvRes <- GetSurvEval(pData = pData, myCum = myHvec, myB = Tres$new_beta, evalT = twopoint)
        # SurvRes_true <- GetSurvEval(pData = pData, myCum = bhest_true, myB = trueBeta, evalT = twopoint)
        # surv1 <- approx(x = TransSurvRes$newCum$T, y = TransSurvRes$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        # survTrue <- approx(x = SurvRes_true$newCum$T, y = SurvRes_true$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        
        rec_beta[i,] <- Tres_btsp$new_beta
        rec_twop_cumH[i,] <- t1
        rec_cumH[i,] <- ttt1
        # rec_twop_surv[i,] <- TransSurvRes$RMST
        # rec_survTrans[i,] <- surv1
    }
    close(pb)
    return(list(rec_beta = rec_beta,
                beta_SE_btsp = apply(rec_beta, 2, sd),
                rec_twop_cumH = rec_twop_cumH,
                cumH_twop_SE_btsp = apply(rec_twop_cumH, 2, sd),
                rec_cumH = rec_cumH,
                cumH_SE_btsp = apply(rec_cumH, 2, sd)))
}

# rec_twop_surv = rec_twop_surv,
# surv_twop_SE_btsp = apply(rec_twop_surv, 2, sd),
# rec_survTrans = rec_survTrans,
# surv_SE_btsp = apply(rec_survTrans, 2, sd))


GetCoverage <- function(estimate, trueBeta,
                        SE) {
    LL <- estimate - 1.96*SE
    UU <- estimate + 1.96*SE
    return(as.numeric(trueBeta<UU & trueBeta>LL))
}


GetSurvEval_RD <- function(pData, myCum, myB, evalT = c(1,2), medianCov = c(50,1,1,0,1,1)) {
    # medianCov = c(median(pData$X1), as.numeric(names(which.max(table(pData$X2)))))
    # medianCov = c(50,1,1,0,1,1)
    newCum = data.frame(H = myCum$hazard * c(exp(medianCov %*% myB)),
                        T = myCum$time)
    newCum$SurvF = exp(-newCum$H)
    newCum <- newCum[order(newCum$T), ]
    newCum <- rbind(c(0,0,1), newCum)
    if(max(newCum$T)<2) {
        newCum <- rbind(newCum, c(max(newCum$H),2,0))
    }
    
    ## get RMST
    func = approxfun(newCum$T, newCum$SurvF)
    RMST = rep(NA, length(evalT))
    for(i in 1:length(evalT)) {
        RMST[i] = integrate(func, 0, evalT[i], subdivisions=2000, stop.on.error = FALSE)$value
    }
    return(list(newCum = newCum,
                RMST = RMST))
}

runBtsp_transCox_RD <- function(aData,
                                pData,
                                twopoint,
                                nbootstrap = 100,
                                trueBeta = c(-0.5, 0.5),
                                cov = c("X1", "X2")) {
    
    rec_beta = rec_twop_cumH = rec_twop_surv <- matrix(NA, nbootstrap, length(cov))
    rec_cumH <- matrix(NA, nbootstrap, 50)
    rec_survTrans = matrix(NA, nbootstrap, 50)
    
    pb = txtProgressBar(min = 0, max = nbootstrap, style = 2, initial = 0) 
    for(i in 1:nbootstrap) {
        setTxtProgressBar(pb, i)
        
        aData_btsp <- aData[sample(1:nrow(aData), nrow(aData), replace = TRUE), ]
        pData_btsp <- pData[sample(1:nrow(pData), nrow(pData), replace = TRUE), ]
        Cout_btsp <- GetAuxSurv(aData_btsp, cov = cov)
        Pout_btsp <- GetPrimaryParam(pData_btsp, q = Cout_btsp$q, estR = Cout_btsp$estR)
        Tres_btsp <- runTransCox_one(Pout_btsp, 
                                     l1 = BestLam1, 
                                     l2 = BestLam2, 
                                     learning_rate = learning_rate, 
                                     nsteps = nsteps,
                                     cov = cov)
        Htime <- Tres_btsp$time
        DupN <- sum(duplicated(Htime)) 
        myHvec <- data.frame(hazard = dQtocumQ(Tres_btsp$new_IntH), time = Htime)
        point50 = seq(1, 120, length.out = 50)
        t1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = twopoint)$y
        ttt1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = point50)$y
        
        # bhest_true = data.frame(time = point50, 
        #                         hazards = point50^2)
        # TransSurvRes <- GetSurvEval_RD(pData = pData, myCum = myHvec, myB = Tres$new_beta, evalT = twopoint)
        # SurvRes_true <- GetSurvEval_RD(pData = pData, myCum = bhest_true, myB = trueBeta, evalT = twopoint)
        # surv1 <- approx(x = TransSurvRes$newCum$T, y = TransSurvRes$newCum$SurvF, xout = point50)$y
        # survTrue <- approx(x = SurvRes_true$newCum$T, y = SurvRes_true$newCum$SurvF, xout = point50)$y
        
        rec_beta[i,] <- Tres_btsp$new_beta
        rec_twop_cumH[i,] <- t1
        rec_cumH[i,] <- ttt1
        # rec_twop_surv[i,] <- TransSurvRes$RMST
        # rec_survTrans[i,] <- surv1
    }
    close(pb)
    return(list(rec_beta = rec_beta,
                beta_SE_btsp = apply(rec_beta, 2, sd),
                rec_twop_cumH = rec_twop_cumH,
                cumH_twop_SE_btsp = apply(rec_twop_cumH, 2, sd),
                rec_cumH = rec_cumH,
                cumH_SE_btsp = apply(rec_cumH, 2, sd)))
}

# medianCov = as.numeric(pData[1,1:2])
GetSurvEval_Cov <- function(pData, myCum, myB, evalT = c(1,2), medianCov = c(0.5,1)) {
    newCum = data.frame(H = myCum$hazard * c(exp(medianCov %*% myB)),
                        T = myCum$time)
    newCum$SurvF = exp(-newCum$H)
    newCum <- newCum[order(newCum$T), ]
    SurvProp <- approx(x = newCum$T, y = newCum$SurvF, xout = twopoint)$y
    return(SurvProp)
}
GetSurvEval_Cov_vec <- function(pData, myCum, myB, evalT = c(1,2)) {
    outProp <- matrix(NA, nrow(pData), length(evalT))
    for(i in 1:nrow(pData)) {
        outProp[i,] <- GetSurvEval_Cov(pData = pData, 
                                       myCum = myCum, 
                                       myB = myB, 
                                       evalT = evalT, 
                                       medianCov = as.numeric(pData[i,1:length(myB)]))
    }
    return(outProp)
}


# rec_twop_surv = rec_twop_surv,
# surv_twop_SE_btsp = apply(rec_twop_surv, 2, sd),
# rec_survTrans = rec_survTrans,
# surv_SE_btsp = apply(rec_survTrans, 2, sd)

