setwd("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")

#### run 100 simulation
#### Setting: main simulation (5 covariates, Unif-Beta in X3)
source("GetData.R")
source("TransCoxRFunctions_auxOnly.R")
library(reticulate)
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))
library(TransCox)
library(survival)
set.seed(1345)

allsetRes <- matrix(NA, 4, 58)
Niter = 100
trueBeta = c(-0.5, 0.5, 0.2, 0.1, 0.1)
biasmat <- matrix(NA, Niter, 12)
Hbiasmat <- matrix(NA, Niter, 8)
RMSTbiasmat <- matrix(NA, Niter, 10)
Hmat <- matrix(NA, Niter, 4)
HdiffMat1 = HdiffMat2 = HdiffMat3 = HdiffMat4 = matrix(NA, Niter, 50)
SurvMat1 = SurvMat2 = SurvMat3 = SurvMat4 = SurvMatTrue = matrix(NA, Niter, 50)
HtwoMat1 = HtwoMat2 = HtwoMat3 = HtwoMat4 = matrix(NA, Niter, 2)

PersPredE1 = PersPredE2 = PersPredE3 = PersPredE4 = matrix(NA, Niter, 2)
AbsPersPredE1 = AbsPersPredE2 = AbsPersPredE3 = AbsPersPredE4 = matrix(NA, Niter, 2)

CumHaz_Trans = CumHaz_Ponly = CumHaz_Both = CumHaz_Str = list()
SurvC_Trans = SurvC_Ponly = SurvC_Both = SurvC_Str = SurvC_True = list()

fullBetaBias = fullHtwoBias = list()
fullPersPred = fullAbsPersPred = list()

nprim = 250
naux = 6400

## The commented part is how we select the learning rate and nstep parameter
# onedata <- GenData_FiveCov(nprim = nprim,
#                            naux = naux,
#                            setting = 1)
# pData = onedata$primData
# aData = onedata$auxData
# LRout <- SelLR_By_BIC(pData, aData, cov = c("X1", "X2", "X3", "X4", "X5"),
#                       lambda1 = 0.1,
#                       lambda2 = 0.1,
#                       learning_rate_vec  = c(0.0005, 0.0008, 0.001, 0.0011, 0.0012, 0.0015, 0.002,0.003),
#                       nsteps_vec = c(700,800,900,1000,1100))

learning_rate = 0.0006  
nsteps = 800
for(nset in 1:4) {
    cat("Setting = ", nset, "\n")
    for(i in 1:Niter) {
        set.seed(i)
        cat(i, ",")
        onedata <- GenData_FiveCov(nprim = nprim,
                                   naux = naux, 
                                   setting = nset,
                                   Dist = "Norm")
        pData = onedata$primData
        aData = onedata$auxData
        
        ###### first, select tuning parameter
        if(i == 1) {
            res <- SelParam_By_BIC(pData, aData, cov = c("X1", "X2", "X3", "X4", "X5"),
                                   statusvar = "status",
                                   lambda1_vec = c(0.1, 0.5, 1, 3, 5, 7),
                                   lambda2_vec = c(0.1, 0.5, 1, 1.5, 2, 3),
                                   learning_rate = learning_rate,
                                   nsteps = nsteps)
            BestLam1 <- res$best_la1
            BestLam2 <- res$best_la2
            cat(BestLam1, ";", BestLam2, "\n")
        }
        
        Cout <- GetAuxSurv(aData, cov = c("X1", "X2", "X3", "X4", "X5"))
        Pout <- GetPrimaryParam(pData, q = Cout$q, estR = Cout$estR)
        
        twopoint = c(0.6, 1.2)
        
        Tres <- runTransCox_one(Pout, 
                                l1 = BestLam1, l2 = BestLam2, learning_rate = learning_rate, 
                                nsteps = nsteps, cov = c("X1", "X2", "X3", "X4", "X5"))
        myHvec <- data.frame(hazard = dQtocumQ(Tres$new_IntH), time = Tres$time)
        TransSurvRes <- GetSurvEval_5Cov(pData = pData, myCum = myHvec, myB = Tres$new_beta, evalT = twopoint)
        T_SurvEst <- GetSurvEval_Cov_vec(pData = pData, myCum = myHvec, myB = Tres$new_beta, evalT = twopoint)
        
        res_Ponly <- coxph(Surv(time, status == 2) ~ X1 + X2 + X3 + X4 + X5, data = pData)
        bhest_Ponly <- basehaz(res_Ponly, centered=FALSE) ## get baseline cumulative hazards
        estR_Ponly <- res_Ponly$coefficients
        SurvRes_Ponly <- GetSurvEval_5Cov(pData = pData, myCum = bhest_Ponly, myB = estR_Ponly, evalT = twopoint)
        Ponly_SurvEst <- GetSurvEval_Cov_vec(pData = pData,myCum = bhest_Ponly, myB = estR_Ponly, evalT = twopoint)
        
        res_both <- coxph(Surv(time, status == 2) ~ X1 + X2 + X3 + X4 + X5, data = rbind(pData, aData))
        bhest_both <- basehaz(res_both, centered=FALSE) ## get baseline cumulative hazards
        estR_both <- res_both$coefficients
        SurvRes_both <- GetSurvEval_5Cov(pData = pData, myCum = bhest_both, myB = estR_both, evalT = twopoint)
        both_SurvEst <- GetSurvEval_Cov_vec(pData = pData,myCum = bhest_both, myB = estR_both, evalT = twopoint)
        
        twodata <- rbind(pData, aData)
        twodata$group <- c(rep(1, nrow(pData)), rep(2, nrow(aData)))
        res_str <- coxph(Surv(time, status == 2) ~ X1 + X2 + X3 + X4 + X5 + strata(group), data = twodata)
        bhest_str <- basehaz(res_str, centered=FALSE) ## get baseline cumulative hazards
        bhest_str <- bhest_str[bhest_str$strata == "group=1", c(1:2)]
        estR_str <- res_str$coefficients
        SurvRes_str <- GetSurvEval_5Cov(pData = pData, myCum = bhest_str, myB = estR_str, evalT = twopoint)
        str_SurvEst <- GetSurvEval_Cov_vec(pData = pData, myCum = bhest_str, myB = estR_str, evalT = twopoint)
        
        point50 = seq(0.15, 2, length.out = 50)^2
        
        bhest_true = data.frame(time = point50, 
                                hazards = point50^2)
        SurvRes_true <- GetSurvEval_5Cov(pData = pData, myCum = bhest_true, myB = trueBeta, evalT = twopoint)
        true_SurvEst <- GetSurvEval_Cov_vec(pData = pData, myCum = bhest_true, myB = trueBeta, evalT = twopoint)
        
        PersPredE1[i,] <- colMeans(T_SurvEst - true_SurvEst)
        PersPredE2[i,] <- colMeans(Ponly_SurvEst - true_SurvEst)
        PersPredE3[i,] <- colMeans(both_SurvEst - true_SurvEst)
        PersPredE4[i,] <- colMeans(str_SurvEst - true_SurvEst)
        
        AbsPersPredE1[i,] <- colMeans(abs(T_SurvEst - true_SurvEst))
        AbsPersPredE2[i,] <- colMeans(abs(Ponly_SurvEst - true_SurvEst))
        AbsPersPredE3[i,] <- colMeans(abs(both_SurvEst - true_SurvEst))
        AbsPersPredE4[i,] <- colMeans(abs(str_SurvEst - true_SurvEst))
        
        t1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = twopoint)$y
        t2 <- approx(x = bhest_Ponly$time, y = bhest_Ponly$hazard, xout = twopoint)$y
        t3 <- approx(x = bhest_both$time, y = bhest_both$hazard, xout = twopoint)$y
        t4 <- approx(x = bhest_str$time, y = bhest_str$hazard, xout = twopoint)$y
        
        ttt1 <- approx(x = myHvec$time, y = myHvec$hazard, xout = seq(0.15, 2, length.out = 50))$y
        ttt2 <- approx(x = bhest_Ponly$time, y = bhest_Ponly$hazard, xout = seq(0.15, 2, length.out = 50))$y
        ttt3 <- approx(x = bhest_both$time, y = bhest_both$hazard, xout = seq(0.15, 2, length.out = 50))$y
        ttt4 <- approx(x = bhest_str$time, y = bhest_str$hazard, xout = seq(0.15, 2, length.out = 50))$y
        
        surv1 <- approx(x = TransSurvRes$newCum$T, y = TransSurvRes$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        surv2 <- approx(x = SurvRes_Ponly$newCum$T, y = SurvRes_Ponly$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        surv3 <- approx(x = SurvRes_both$newCum$T, y = SurvRes_both$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        surv4 <- approx(x = SurvRes_str$newCum$T, y = SurvRes_str$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        survTrue <- approx(x = SurvRes_true$newCum$T, y = SurvRes_true$newCum$SurvF, xout = seq(0, 2, length.out = 50))$y
        
        threeBias <- c(Tres$new_beta[1:3] - trueBeta[1:3],
                       estR_Ponly[1:3] - trueBeta[1:3],
                       estR_both[1:3] - trueBeta[1:3],
                       estR_str[1:3] - trueBeta[1:3])
        
        HtwoMat1[i, ] <- t1 - twopoint^2
        HtwoMat2[i, ] <- t2 - twopoint^2
        HtwoMat3[i, ] <- t3 - twopoint^2
        HtwoMat4[i, ] <- t4 - twopoint^2
        
        HdiffMat1[i, ] <- ttt1
        HdiffMat2[i, ] <- ttt2
        HdiffMat3[i, ] <- ttt3
        HdiffMat4[i, ] <- ttt4
        
        SurvMat1[i, ] <- surv1
        SurvMat2[i, ] <- surv2
        SurvMat3[i, ] <- surv3
        SurvMat4[i, ] <- surv4
        
        SurvMatTrue[i, ] <- survTrue
        
        RMSTBias <- c(TransSurvRes$RMST - SurvRes_true$RMST,
                      SurvRes_Ponly$RMST - SurvRes_true$RMST,
                      SurvRes_both$RMST - SurvRes_true$RMST,
                      SurvRes_str$RMST - SurvRes_true$RMST,
                      SurvRes_true$RMST)
        
        biasmat[i,] <- threeBias
        RMSTbiasmat[i, ] <- RMSTBias
        Hbiasmat[i,] <- c(t1 - twopoint^2,
                          t2 - twopoint^2,
                          t3 - twopoint^2,
                          t4 - twopoint^2)
    }
    
    mean_vec = signif(colMeans(biasmat), 4)
    sd_vec = signif(apply(biasmat, 2, sd), 4)
    Hmat_mean = signif(c(colMeans(HtwoMat1, na.rm = TRUE),
                         colMeans(HtwoMat2, na.rm = TRUE),
                         colMeans(HtwoMat3, na.rm = TRUE),
                         colMeans(HtwoMat4, na.rm = TRUE)), 4)
    Hmat_sd = signif(c(apply(HtwoMat1, 2, sd, na.rm = TRUE),
                       apply(HtwoMat2, 2, sd, na.rm = TRUE),
                       apply(HtwoMat3, 2, sd, na.rm = TRUE),
                       apply(HtwoMat4, 2, sd, na.rm = TRUE)), 4)
    RMST_mean = signif(colMeans(RMSTbiasmat), 4)
    RMST_sd = signif(apply(RMSTbiasmat, 2, sd), 4)
    
    PersPredMat <- cbind(PersPredE1, PersPredE2, PersPredE3, PersPredE4)
    PP_mean <- signif(colMeans(PersPredMat), 4)
    PP_sd <- signif(apply(PersPredMat, 2, sd), 4)
    AbsPersPredMat <- cbind(AbsPersPredE1, AbsPersPredE2, AbsPersPredE3, AbsPersPredE4)
    APP_mean <- signif(colMeans(AbsPersPredMat), 4)
    APP_sd <- signif(apply(AbsPersPredMat, 2, sd), 4)
    
    allsetRes[nset, ] <- c(paste0(mean_vec, "(", sd_vec, ")"),
                           signif(colMeans(biasmat)^2+apply(biasmat, 2, var), 4),
                           c(paste0(Hmat_mean, "(", Hmat_sd, ")")),
                           c(paste0(RMST_mean, "(", RMST_sd, ")")),
                           c(paste0(PP_mean, "(", PP_sd, ")")),
                           c(paste0(APP_mean, "(", APP_sd, ")")))
    
    fullPersPred[[nset]] = PersPredMat
    fullAbsPersPred[[nset]] = AbsPersPredMat
    
    CumHaz_Trans[[nset]] <- HdiffMat1
    CumHaz_Ponly[[nset]] <- HdiffMat2
    CumHaz_Both[[nset]] <- HdiffMat3
    CumHaz_Str[[nset]] <- HdiffMat4
    
    SurvC_Trans[[nset]] <- SurvMat1
    SurvC_Ponly[[nset]] <- SurvMat2
    SurvC_Both[[nset]] <- SurvMat3
    SurvC_Str[[nset]] <- SurvMat4
    
    SurvC_True[[nset]] <- SurvMatTrue
    
    fullBetaBias[[nset]] <- biasmat
    fullHtwoBias[[nset]] <- Hbiasmat
}


colnames(allsetRes) <- c("TransCox_b1_bias (sd)", "TransCox_b2_bias (sd)", "TransCox_b3_bias (sd)",
                         "CoxP_b1_bias (sd)", "CoxP_b2_bias (sd)", "CoxP_b3_bias (sd)",
                         "CoxB_b1_bias", "CoxB_b2_bias", "CoxB_b3_bias (sd)",
                         "CoxStr_b1_bias", "CoxStr_b2_bias", "CoxStr_b3_bias (sd)",
                         "TransCox_b1_MSE", "TransCox_b2_MSE", "TransCox_b3_MSE",
                         "CoxP_b1_MSE", "CoxP_b2_MSE", "CoxP_b3_MSE",
                         "CoxB_b1_MSE", "CoxB_b2_MSE", "CoxB_b3_MSE",
                         "CoxStr_b1_MSE", "CoxStr_b2_MSE", "CoxStr_b3_MSE",
                         "TransCox_meanH_bias_1 (sd)", "TransCox_meanH_bias_2 (sd)", 
                         "CoxP_meanH_bias_1 (sd)", "CoxP_meanH_bias_2 (sd)",
                         "CoxB_meanH_bias_1 (sd)", "CoxB_meanH_bias_2 (sd)",
                         "CoxStr_meanH_bias_1 (sd)", "CoxStr_meanH_bias_2 (sd)",
                         "TransCox_RMST_bias_1 (sd)", "TransCox_RMST_bias_2 (sd)", 
                         "CoxP_RMST_bias_1 (sd)", "CoxP_RMST_bias_2 (sd)",
                         "CoxB_RMST_bias_1 (sd)", "CoxB_RMST_bias_2 (sd)",
                         "CoxStr_RMST_bias_1 (sd)", "CoxStr_RMST_bias_2 (sd)",
                         "TRUE_RMST_bias_1 (sd)", "TRUE_RMST_bias_2 (sd)",
                         "TransCox_PP_bias_1 (sd)", "TransCox_PP_bias_2 (sd)", 
                         "CoxP_PP_bias_1 (sd)", "CoxP_P_bias_2 (sd)",
                         "CoxB_PP_bias_1 (sd)", "CoxB_PP_bias_2 (sd)",
                         "CoxStr_PP_bias_1 (sd)", "CoxStr_PP_bias_2 (sd)",
                         "TransCox_APP_bias_1 (sd)", "TransCox_APP_bias_2 (sd)", 
                         "CoxP_APP_bias_1 (sd)", "CoxP_APP_bias_2 (sd)",
                         "CoxB_APP_bias_1 (sd)", "CoxB_APP_bias_2 (sd)",
                         "CoxStr_APP_bias_1 (sd)", "CoxStr_APP_bias_2 (sd)")

rownames(allsetRes) <- paste0("Setting ", 1:4)

allsetRes

save(allsetRes, CumHaz_Trans,
     CumHaz_Ponly, CumHaz_Both,
     SurvC_Trans, SurvC_Ponly,
     SurvC_Both, SurvC_True,
     fullBetaBias, fullHtwoBias,
     fullPersPred, fullAbsPersPred,
     file = paste0("Cov5_Normal_SimResults_", nprim, "_", naux, ".RData"))

